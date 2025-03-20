// word width define
`define WI 9
`define WF 26
`define WFPHASE 16
`define WCNT 10

// -------------------------------------------------------
// Module Name: LO_PSYNC_CTRL
// Function: measure lo i/q phase according to 2bit sample data, generate 1bit signal to exchange lo i/q phase
// Author: Yang Yumeng Date: 4/12 2023
// Version: v1p0
// -------------------------------------------------------
module LOdiv2_PSYNC_CTRL (
NARST,
EN,
CLK,
SYS_REF,
PLL_FCW,
LO_DIV,
LO_STATE,
DONE,
EXCHOP
);

// io decliration
input NARST; // asynchronouse nrst
input EN; 
input CLK; // trigger at posedge
input SYS_REF; // system ref
input [9+26-1:0] PLL_FCW; // same as the FCW sent to PLL
input [2:0] LO_DIV; // div 2/4/6/8/16/32/64/128/256
input [1:0] LO_STATE; // sampled lo phase
output reg [1:0] DONE; // phase synchronize fsm state sign; [1] 0: execute, 1: done; [0] 0: failed, 1: succeed
output reg EXCHOP; // exchange i/o operation signal, 0 needn't exchange, 1 exchange

// fsm param define
localparam fsm_state_measure = 2'b00; // measure & calculate phase error
localparam fsm_state_genexch = 2'b01; // generate EXCHOP & validate
localparam fsm_state_done = 2'b10; // done

// internal signal
wire NRST;
reg sys_ref_d1;
reg sys_ref_d2;
reg sys_ref_d3;
reg sys_ref_d4;
wire sys_comb;
wire sys_ctrl;
reg sys_psync_en;
wire [9+26-1:0] LO_FCW_F;
wire [9+26-1:0] PACCUM_LIMIT; // LO digital freq accumulation limitation.
reg PACCUM_c;
reg [9+26-1:0] PACCUM_s;
wire [3:0] LO_DIV_add;
wire [26PHASE-1:0] dphase_lo_c;
reg [26PHASE-1:0] dphase_lo_m;
reg [10-1:0] window;
reg flag;
wire [26PHASE-1:0] diffphase;// (ufix), 0+16, [0deg,360deg)
reg [26PHASE-1+10:0] diffphase_sum;
wire [26PHASE-1:0] diffphase_avg;
wire [26PHASE-1:0] diffphase_shift;
wire vote;
reg [1:0] fsm_state;

// code begin

SYNCRSTGEN U0_SYNCRST( .CLK (CLK), .NARST (NARST), .NRST (NRST), .NRST1 (), .NRST2());

// generate phase synchronouse enable signal
always @ (posedge CLK) begin
	sys_ref_d1 <= SYS_REF;
	sys_ref_d2 <= sys_ref_d1;
	sys_ref_d3 <= sys_ref_d2;
	sys_ref_d4 <= sys_ref_d3;
end

assign sys_comb = sys_ref_d3 & (~sys_ref_d4);
assign sys_ctrl = EN & sys_comb;

always @ (posedge CLK or negedge NRST) begin
	if (!NRST) begin 
		sys_psync_en <= 1'b0;
	end else begin
		if ((sys_psync_en==1'b0)&&(sys_ctrl==1'b1)) begin
			sys_psync_en <= 1'b1;
		end else if (EN==1'b0) begin
			sys_psync_en <= 1'b0;
		end
	end
end

// digital phase generator
// lo freq calc
assign LO_DIV_add = LO_DIV;
assign LO_FCW_F = PLL_FCW - ((PLL_FCW >> (26+LO_DIV_add)) << (26+LO_DIV_add));
assign PACCUM_LIMIT = 1'b1 << (26+LO_DIV_add);
assign dphase_lo_c = (|((1'b1<<(26-26PHASE-1+LO_DIV_add))&PACCUM_s))? ((PACCUM_s>>(26-26PHASE+LO_DIV_add)) + 1'b1): ((PACCUM_s>>(26-26PHASE+LO_DIV_add)) + 1'b0);

always @ (posedge CLK or negedge NRST) begin
	if (!NRST) begin
		PACCUM_s <= 0;
	end else if (sys_psync_en) begin
		if ( ( PACCUM_s + LO_FCW_F ) < PACCUM_LIMIT ) begin
			PACCUM_s <= PACCUM_s + LO_FCW_F;
		end else begin
			PACCUM_s <= PACCUM_s + LO_FCW_F - PACCUM_LIMIT;
		end
	end else begin
		PACCUM_s <= 0;
	end
end

// map lo state to quantization analog phase
// sampler
always @ (posedge CLK or negedge NRST) begin
	if (!NRST) dphase_lo_m <= 0;
	else if (sys_psync_en) begin
		case (LO_STATE)
			2'b10: dphase_lo_m <= 0;
			2'b11: dphase_lo_m <= {2'b01, {(26PHASE-2){1'b0}}}; // 090 deg
			2'b01: dphase_lo_m <= {2'b10, {(26PHASE-2){1'b0}}}; // 180 deg
			2'b00: dphase_lo_m <= {2'b11, {(26PHASE-2){1'b0}}}; // 270 deg
		endcase
	end else begin
		dphase_lo_m <= 0;
	end
end

assign diffphase = dphase_lo_m - dphase_lo_c;
assign diffphase_avg = diffphase_sum >> (10);

always @ (posedge CLK or negedge NRST) begin
	if (!NRST) begin
		diffphase_sum <= 0;
	end else if (sys_psync_en) begin
		diffphase_sum <= flag? 0: (diffphase_sum + {{(10){diffphase[26PHASE-1]}}, diffphase});
	end
end

// vote 
assign diffphase_shift = diffphase_avg - {2'b01, {(26PHASE-2){1'b0}}};
assign vote = diffphase_shift[26PHASE-1];

// phase error calculation

// generate 2048 cycles window
always @ (posedge CLK or negedge NRST) begin
	if (!NRST) begin
		{flag, window} <= 0;
	end else if (sys_psync_en) begin
		{flag, window} <= window + 1'b1;
	end
end

// state machine
// generate exchange i/o operation signal

always @ (posedge CLK or negedge NRST) begin
	if (!NRST) begin
		fsm_state <= 0;
		DONE <= 0;
		EXCHOP <= 0;
	end else if (sys_psync_en&flag) begin
		case (fsm_state)
			fsm_state_measure: begin
				fsm_state <= fsm_state_genexch;
				DONE <= 2'b00;
				EXCHOP <= vote;
			end
			fsm_state_genexch: begin
				fsm_state <= fsm_state_done;
				DONE <= {1'b1, ~vote};
				EXCHOP <= EXCHOP;
			end
			fsm_state_done: begin
				fsm_state <= fsm_state_done;
				DONE <= DONE;
				EXCHOP <= EXCHOP;
			end
		endcase
	end
end

// // test
// integer fp1;

// initial begin
	// fp1 = $fopen("./diffphase.txt");
// end

// always @ (posedge flag) begin
	// $fstrobe(fp1, "%d", $unsigned(diffphase_avg));
// end

endmodule