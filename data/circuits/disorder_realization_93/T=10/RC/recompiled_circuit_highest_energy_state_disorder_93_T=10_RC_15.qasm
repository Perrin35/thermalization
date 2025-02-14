OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2095069) q[0];
sx q[0];
rz(-1.0340438) q[0];
sx q[0];
rz(0.75420585) q[0];
rz(1.7475313) q[1];
sx q[1];
rz(-0.59352195) q[1];
sx q[1];
rz(-1.4374179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2090354) q[0];
sx q[0];
rz(-0.13929312) q[0];
sx q[0];
rz(-2.7053231) q[0];
rz(1.1815156) q[2];
sx q[2];
rz(-1.6061651) q[2];
sx q[2];
rz(0.25098342) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9086409) q[1];
sx q[1];
rz(-0.65083069) q[1];
sx q[1];
rz(-2.6292168) q[1];
rz(-3.0079477) q[3];
sx q[3];
rz(-1.9507532) q[3];
sx q[3];
rz(-2.8524672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6393911) q[2];
sx q[2];
rz(-2.907967) q[2];
sx q[2];
rz(2.8761253) q[2];
rz(-1.4207077) q[3];
sx q[3];
rz(-1.1666433) q[3];
sx q[3];
rz(1.407912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0625286) q[0];
sx q[0];
rz(-1.205235) q[0];
sx q[0];
rz(-2.089654) q[0];
rz(-2.1417292) q[1];
sx q[1];
rz(-1.544416) q[1];
sx q[1];
rz(2.3725407) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3628415) q[0];
sx q[0];
rz(-1.0386416) q[0];
sx q[0];
rz(2.6274985) q[0];
rz(-pi) q[1];
rz(-3.0508383) q[2];
sx q[2];
rz(-2.5620775) q[2];
sx q[2];
rz(-1.3578292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96892541) q[1];
sx q[1];
rz(-1.7255867) q[1];
sx q[1];
rz(2.0230146) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6399986) q[3];
sx q[3];
rz(-2.6644649) q[3];
sx q[3];
rz(-2.9692544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.179788) q[2];
sx q[2];
rz(-2.3369117) q[2];
sx q[2];
rz(1.5847607) q[2];
rz(0.70476091) q[3];
sx q[3];
rz(-2.9289991) q[3];
sx q[3];
rz(1.0889277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3253118) q[0];
sx q[0];
rz(-0.75958696) q[0];
sx q[0];
rz(-2.3749206) q[0];
rz(0.39241544) q[1];
sx q[1];
rz(-2.2800443) q[1];
sx q[1];
rz(-0.79889417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5094389) q[0];
sx q[0];
rz(-1.5804133) q[0];
sx q[0];
rz(-2.0548178) q[0];
x q[1];
rz(0.33943265) q[2];
sx q[2];
rz(-2.8763874) q[2];
sx q[2];
rz(-2.3496534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0943999) q[1];
sx q[1];
rz(-0.59941888) q[1];
sx q[1];
rz(1.0034189) q[1];
x q[2];
rz(1.6774936) q[3];
sx q[3];
rz(-1.7472526) q[3];
sx q[3];
rz(0.7430232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7400292) q[2];
sx q[2];
rz(-0.26569685) q[2];
sx q[2];
rz(-0.024988739) q[2];
rz(1.7449024) q[3];
sx q[3];
rz(-1.9091505) q[3];
sx q[3];
rz(-0.11740824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5806737) q[0];
sx q[0];
rz(-0.54045254) q[0];
sx q[0];
rz(-0.29676357) q[0];
rz(-0.98776039) q[1];
sx q[1];
rz(-1.2514273) q[1];
sx q[1];
rz(-2.3938649) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0714617) q[0];
sx q[0];
rz(-0.80118766) q[0];
sx q[0];
rz(1.8772932) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1224611) q[2];
sx q[2];
rz(-0.90662642) q[2];
sx q[2];
rz(-2.3519206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8104483) q[1];
sx q[1];
rz(-0.40256009) q[1];
sx q[1];
rz(3.1123575) q[1];
rz(-pi) q[2];
rz(-2.4676855) q[3];
sx q[3];
rz(-2.0567937) q[3];
sx q[3];
rz(0.35714071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.957754) q[2];
sx q[2];
rz(-2.2695393) q[2];
sx q[2];
rz(-2.6217065) q[2];
rz(-2.1947491) q[3];
sx q[3];
rz(-1.1917453) q[3];
sx q[3];
rz(0.051518353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87676048) q[0];
sx q[0];
rz(-0.89986372) q[0];
sx q[0];
rz(-1.1967891) q[0];
rz(-1.4747249) q[1];
sx q[1];
rz(-0.80497545) q[1];
sx q[1];
rz(-2.8428452) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099649053) q[0];
sx q[0];
rz(-1.9199326) q[0];
sx q[0];
rz(-0.70167244) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0709956) q[2];
sx q[2];
rz(-0.83865863) q[2];
sx q[2];
rz(3.0867019) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50089624) q[1];
sx q[1];
rz(-1.4569893) q[1];
sx q[1];
rz(0.76687419) q[1];
rz(-1.1029873) q[3];
sx q[3];
rz(-2.3846979) q[3];
sx q[3];
rz(1.6136957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3195191) q[2];
sx q[2];
rz(-0.75054979) q[2];
sx q[2];
rz(3.1263515) q[2];
rz(-3.0139253) q[3];
sx q[3];
rz(-2.7805507) q[3];
sx q[3];
rz(2.291919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3540038) q[0];
sx q[0];
rz(-2.6241527) q[0];
sx q[0];
rz(-2.2678243) q[0];
rz(1.9173701) q[1];
sx q[1];
rz(-1.0226378) q[1];
sx q[1];
rz(-1.9220985) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50724678) q[0];
sx q[0];
rz(-2.0197045) q[0];
sx q[0];
rz(2.6953791) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72645006) q[2];
sx q[2];
rz(-1.9973576) q[2];
sx q[2];
rz(-0.62131617) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4095597) q[1];
sx q[1];
rz(-0.28342208) q[1];
sx q[1];
rz(0.83614142) q[1];
x q[2];
rz(0.55762724) q[3];
sx q[3];
rz(-2.4708797) q[3];
sx q[3];
rz(-2.022555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.445861) q[2];
sx q[2];
rz(-1.4706688) q[2];
sx q[2];
rz(-2.640558) q[2];
rz(1.3387574) q[3];
sx q[3];
rz(-0.41156358) q[3];
sx q[3];
rz(1.1494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5981136) q[0];
sx q[0];
rz(-1.9444332) q[0];
sx q[0];
rz(0.57394779) q[0];
rz(-2.5698938) q[1];
sx q[1];
rz(-1.0400925) q[1];
sx q[1];
rz(-1.5843102) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1645419) q[0];
sx q[0];
rz(-0.91616154) q[0];
sx q[0];
rz(-0.47619168) q[0];
rz(1.1102067) q[2];
sx q[2];
rz(-0.71515036) q[2];
sx q[2];
rz(0.88727027) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3347335) q[1];
sx q[1];
rz(-2.1424531) q[1];
sx q[1];
rz(2.566659) q[1];
rz(-1.9628223) q[3];
sx q[3];
rz(-1.2268664) q[3];
sx q[3];
rz(1.0052538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.456936) q[2];
sx q[2];
rz(-1.9921649) q[2];
sx q[2];
rz(-1.9880902) q[2];
rz(1.614511) q[3];
sx q[3];
rz(-1.0228415) q[3];
sx q[3];
rz(1.3326741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.033878) q[0];
sx q[0];
rz(-2.3623473) q[0];
sx q[0];
rz(-0.23767924) q[0];
rz(-2.4029845) q[1];
sx q[1];
rz(-1.9269201) q[1];
sx q[1];
rz(0.019066378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4486074) q[0];
sx q[0];
rz(-3.0685022) q[0];
sx q[0];
rz(1.4574217) q[0];
rz(2.4408728) q[2];
sx q[2];
rz(-2.2402175) q[2];
sx q[2];
rz(0.6159516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42603675) q[1];
sx q[1];
rz(-1.6095248) q[1];
sx q[1];
rz(1.6524131) q[1];
rz(-pi) q[2];
rz(-0.98293368) q[3];
sx q[3];
rz(-2.1275828) q[3];
sx q[3];
rz(-2.9012321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0154524) q[2];
sx q[2];
rz(-1.3498053) q[2];
sx q[2];
rz(-2.9665185) q[2];
rz(1.3779047) q[3];
sx q[3];
rz(-0.62052369) q[3];
sx q[3];
rz(-0.11848816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058067583) q[0];
sx q[0];
rz(-1.7836934) q[0];
sx q[0];
rz(2.626626) q[0];
rz(-0.99634755) q[1];
sx q[1];
rz(-1.0641655) q[1];
sx q[1];
rz(-1.6641585) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9148902) q[0];
sx q[0];
rz(-0.88719598) q[0];
sx q[0];
rz(-0.76163389) q[0];
rz(-pi) q[1];
rz(-3.0918432) q[2];
sx q[2];
rz(-0.58895117) q[2];
sx q[2];
rz(2.4204554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.094992735) q[1];
sx q[1];
rz(-2.0741724) q[1];
sx q[1];
rz(-0.51790389) q[1];
rz(-pi) q[2];
rz(-1.1682801) q[3];
sx q[3];
rz(-0.63197836) q[3];
sx q[3];
rz(3.0807854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.065206334) q[2];
sx q[2];
rz(-1.6519158) q[2];
sx q[2];
rz(-0.15929407) q[2];
rz(-1.6431036) q[3];
sx q[3];
rz(-2.4703333) q[3];
sx q[3];
rz(-1.338965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945187) q[0];
sx q[0];
rz(-3.0975332) q[0];
sx q[0];
rz(-2.2799168) q[0];
rz(1.8623976) q[1];
sx q[1];
rz(-2.8586614) q[1];
sx q[1];
rz(-2.9544725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2047855) q[0];
sx q[0];
rz(-1.7996368) q[0];
sx q[0];
rz(2.8605927) q[0];
rz(-2.6473706) q[2];
sx q[2];
rz(-1.3307327) q[2];
sx q[2];
rz(0.46948813) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4315553) q[1];
sx q[1];
rz(-1.3042684) q[1];
sx q[1];
rz(0.73898594) q[1];
x q[2];
rz(0.6660717) q[3];
sx q[3];
rz(-2.6544673) q[3];
sx q[3];
rz(-2.7745326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6992135) q[2];
sx q[2];
rz(-2.1670161) q[2];
sx q[2];
rz(-1.7787735) q[2];
rz(1.6393225) q[3];
sx q[3];
rz(-2.5876741) q[3];
sx q[3];
rz(-1.7207918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61158553) q[0];
sx q[0];
rz(-1.3349608) q[0];
sx q[0];
rz(0.0064749574) q[0];
rz(-2.6352128) q[1];
sx q[1];
rz(-0.19229278) q[1];
sx q[1];
rz(-2.2813588) q[1];
rz(-1.8851595) q[2];
sx q[2];
rz(-0.4421352) q[2];
sx q[2];
rz(3.0015702) q[2];
rz(-0.68001775) q[3];
sx q[3];
rz(-2.802805) q[3];
sx q[3];
rz(-1.203385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
