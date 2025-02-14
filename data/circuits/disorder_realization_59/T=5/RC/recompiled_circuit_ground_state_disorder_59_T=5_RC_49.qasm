OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6049603) q[0];
sx q[0];
rz(-2.9807615) q[0];
sx q[0];
rz(0.39128006) q[0];
rz(-0.04303509) q[1];
sx q[1];
rz(4.7624762) q[1];
sx q[1];
rz(12.227992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94253263) q[0];
sx q[0];
rz(-1.3056268) q[0];
sx q[0];
rz(1.2808179) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0909371) q[2];
sx q[2];
rz(-1.7850375) q[2];
sx q[2];
rz(3.0681477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8336645) q[1];
sx q[1];
rz(-0.90667716) q[1];
sx q[1];
rz(-0.2425027) q[1];
rz(-pi) q[2];
rz(-1.4621327) q[3];
sx q[3];
rz(-1.9851284) q[3];
sx q[3];
rz(1.765583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7732064) q[2];
sx q[2];
rz(-0.97969222) q[2];
sx q[2];
rz(-1.0580019) q[2];
rz(-3.0392785) q[3];
sx q[3];
rz(-1.5240069) q[3];
sx q[3];
rz(-1.4003096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8232987) q[0];
sx q[0];
rz(-0.75495356) q[0];
sx q[0];
rz(0.21389432) q[0];
rz(-1.8077883) q[1];
sx q[1];
rz(-2.6413481) q[1];
sx q[1];
rz(2.852829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28604668) q[0];
sx q[0];
rz(-0.65546765) q[0];
sx q[0];
rz(1.2857057) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83728055) q[2];
sx q[2];
rz(-1.0425261) q[2];
sx q[2];
rz(2.3579896) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3067291) q[1];
sx q[1];
rz(-2.9188061) q[1];
sx q[1];
rz(1.7053717) q[1];
x q[2];
rz(0.43536868) q[3];
sx q[3];
rz(-1.8791128) q[3];
sx q[3];
rz(-1.8086019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5220962) q[2];
sx q[2];
rz(-2.3213826) q[2];
sx q[2];
rz(1.8572469) q[2];
rz(1.8340825) q[3];
sx q[3];
rz(-1.8239572) q[3];
sx q[3];
rz(-0.6984843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73806015) q[0];
sx q[0];
rz(-2.0750676) q[0];
sx q[0];
rz(2.2962978) q[0];
rz(2.2360133) q[1];
sx q[1];
rz(-0.60494345) q[1];
sx q[1];
rz(-1.9639429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15525165) q[0];
sx q[0];
rz(-3.0444859) q[0];
sx q[0];
rz(-2.9349021) q[0];
x q[1];
rz(-2.6180116) q[2];
sx q[2];
rz(-1.0835907) q[2];
sx q[2];
rz(-0.52766818) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78598394) q[1];
sx q[1];
rz(-2.5289383) q[1];
sx q[1];
rz(-0.36566465) q[1];
rz(-pi) q[2];
rz(-1.3602474) q[3];
sx q[3];
rz(-0.65587924) q[3];
sx q[3];
rz(-2.0149503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5251069) q[2];
sx q[2];
rz(-2.1973886) q[2];
sx q[2];
rz(0.13060972) q[2];
rz(1.3579926) q[3];
sx q[3];
rz(-1.0087174) q[3];
sx q[3];
rz(-1.85359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4865049) q[0];
sx q[0];
rz(-0.26432744) q[0];
sx q[0];
rz(-0.41410145) q[0];
rz(-0.86530238) q[1];
sx q[1];
rz(-1.9046141) q[1];
sx q[1];
rz(0.6032595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7287528) q[0];
sx q[0];
rz(-1.1825292) q[0];
sx q[0];
rz(-0.30465841) q[0];
x q[1];
rz(-2.7895945) q[2];
sx q[2];
rz(-1.504045) q[2];
sx q[2];
rz(2.091832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8131667) q[1];
sx q[1];
rz(-0.35086497) q[1];
sx q[1];
rz(3.0163832) q[1];
rz(-pi) q[2];
rz(-0.11988414) q[3];
sx q[3];
rz(-1.2601687) q[3];
sx q[3];
rz(1.436123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6174751) q[2];
sx q[2];
rz(-1.6092499) q[2];
sx q[2];
rz(-2.4821846) q[2];
rz(0.13255969) q[3];
sx q[3];
rz(-2.5949251) q[3];
sx q[3];
rz(0.70422712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7947614) q[0];
sx q[0];
rz(-2.1857388) q[0];
sx q[0];
rz(-2.8723248) q[0];
rz(1.6412093) q[1];
sx q[1];
rz(-1.2790479) q[1];
sx q[1];
rz(1.0795116) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3985232) q[0];
sx q[0];
rz(-2.4691104) q[0];
sx q[0];
rz(0.29015707) q[0];
x q[1];
rz(2.4731878) q[2];
sx q[2];
rz(-2.2434017) q[2];
sx q[2];
rz(2.2318537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36240807) q[1];
sx q[1];
rz(-1.0917062) q[1];
sx q[1];
rz(-0.37491231) q[1];
x q[2];
rz(0.57860878) q[3];
sx q[3];
rz(-2.7670585) q[3];
sx q[3];
rz(2.7079564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.72011224) q[2];
sx q[2];
rz(-2.5613027) q[2];
sx q[2];
rz(0.87023467) q[2];
rz(-1.8057711) q[3];
sx q[3];
rz(-1.8060874) q[3];
sx q[3];
rz(-0.67572063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0077591) q[0];
sx q[0];
rz(-0.02820153) q[0];
sx q[0];
rz(2.5303685) q[0];
rz(0.73355833) q[1];
sx q[1];
rz(-1.4708142) q[1];
sx q[1];
rz(0.38331097) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0548504) q[0];
sx q[0];
rz(-0.87854687) q[0];
sx q[0];
rz(1.9350497) q[0];
rz(0.82742847) q[2];
sx q[2];
rz(-1.7199923) q[2];
sx q[2];
rz(2.0713446) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94180471) q[1];
sx q[1];
rz(-2.2082303) q[1];
sx q[1];
rz(0.15934039) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72630693) q[3];
sx q[3];
rz(-1.8890813) q[3];
sx q[3];
rz(2.4232211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7395301) q[2];
sx q[2];
rz(-0.5548839) q[2];
sx q[2];
rz(2.128672) q[2];
rz(-0.5976451) q[3];
sx q[3];
rz(-1.7326071) q[3];
sx q[3];
rz(-2.2156318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0585854) q[0];
sx q[0];
rz(-0.42917955) q[0];
sx q[0];
rz(0.035932628) q[0];
rz(1.9195456) q[1];
sx q[1];
rz(-0.67600328) q[1];
sx q[1];
rz(-2.8935208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1897514) q[0];
sx q[0];
rz(-1.763367) q[0];
sx q[0];
rz(1.309881) q[0];
rz(-pi) q[1];
rz(-0.49068292) q[2];
sx q[2];
rz(-1.4457189) q[2];
sx q[2];
rz(-0.016029257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2100468) q[1];
sx q[1];
rz(-1.7341053) q[1];
sx q[1];
rz(-1.1985267) q[1];
rz(-pi) q[2];
rz(-1.5682427) q[3];
sx q[3];
rz(-1.6757343) q[3];
sx q[3];
rz(1.6014639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0849358) q[2];
sx q[2];
rz(-2.0818905) q[2];
sx q[2];
rz(-2.6094931) q[2];
rz(0.45446864) q[3];
sx q[3];
rz(-1.5458509) q[3];
sx q[3];
rz(-1.8475629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249975) q[0];
sx q[0];
rz(-1.0497365) q[0];
sx q[0];
rz(-2.9014034) q[0];
rz(0.60513085) q[1];
sx q[1];
rz(-1.3477707) q[1];
sx q[1];
rz(2.4952707) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511348) q[0];
sx q[0];
rz(-0.86800985) q[0];
sx q[0];
rz(2.7443307) q[0];
rz(-pi) q[1];
rz(-0.036810151) q[2];
sx q[2];
rz(-1.9318077) q[2];
sx q[2];
rz(-1.6474468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5390215) q[1];
sx q[1];
rz(-2.638431) q[1];
sx q[1];
rz(-0.43518592) q[1];
rz(2.2262832) q[3];
sx q[3];
rz(-1.743225) q[3];
sx q[3];
rz(0.61911914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1120844) q[2];
sx q[2];
rz(-1.6479475) q[2];
sx q[2];
rz(-3.0599111) q[2];
rz(-0.84780848) q[3];
sx q[3];
rz(-2.732328) q[3];
sx q[3];
rz(2.9186287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001556) q[0];
sx q[0];
rz(-2.2969022) q[0];
sx q[0];
rz(2.2499114) q[0];
rz(1.2314697) q[1];
sx q[1];
rz(-2.6530118) q[1];
sx q[1];
rz(-2.5317392) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632034) q[0];
sx q[0];
rz(-2.0318306) q[0];
sx q[0];
rz(1.0151063) q[0];
x q[1];
rz(-0.78471176) q[2];
sx q[2];
rz(-1.5605188) q[2];
sx q[2];
rz(-1.0148054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5130723) q[1];
sx q[1];
rz(-1.2781464) q[1];
sx q[1];
rz(2.3696632) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9159007) q[3];
sx q[3];
rz(-1.5865579) q[3];
sx q[3];
rz(1.0712766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2764728) q[2];
sx q[2];
rz(-2.1742994) q[2];
sx q[2];
rz(-2.8709732) q[2];
rz(-0.5451777) q[3];
sx q[3];
rz(-1.3278278) q[3];
sx q[3];
rz(2.8934208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4973064) q[0];
sx q[0];
rz(-1.8338642) q[0];
sx q[0];
rz(0.18648952) q[0];
rz(-1.9322152) q[1];
sx q[1];
rz(-1.3554363) q[1];
sx q[1];
rz(-3.1219416) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.540003) q[0];
sx q[0];
rz(-2.181291) q[0];
sx q[0];
rz(0.50525804) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80060963) q[2];
sx q[2];
rz(-2.7749546) q[2];
sx q[2];
rz(-3.0406102) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33981048) q[1];
sx q[1];
rz(-0.385901) q[1];
sx q[1];
rz(-0.41706155) q[1];
x q[2];
rz(2.2084633) q[3];
sx q[3];
rz(-2.2266042) q[3];
sx q[3];
rz(2.3011014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43202117) q[2];
sx q[2];
rz(-0.95546466) q[2];
sx q[2];
rz(1.552399) q[2];
rz(-3.0324557) q[3];
sx q[3];
rz(-1.8692632) q[3];
sx q[3];
rz(-0.2763589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6266262) q[0];
sx q[0];
rz(-1.4780541) q[0];
sx q[0];
rz(1.9274101) q[0];
rz(-1.7264438) q[1];
sx q[1];
rz(-2.1198004) q[1];
sx q[1];
rz(-1.459495) q[1];
rz(-1.8564275) q[2];
sx q[2];
rz(-2.0002596) q[2];
sx q[2];
rz(3.132435) q[2];
rz(0.23701238) q[3];
sx q[3];
rz(-2.7494299) q[3];
sx q[3];
rz(1.3361479) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
