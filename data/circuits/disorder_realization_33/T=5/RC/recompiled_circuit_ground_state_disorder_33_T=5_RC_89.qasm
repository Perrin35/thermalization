OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6991718) q[0];
sx q[0];
rz(-0.81096634) q[0];
sx q[0];
rz(0.45642689) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(-0.95093095) q[1];
sx q[1];
rz(-0.18263291) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4066494) q[0];
sx q[0];
rz(-1.8514086) q[0];
sx q[0];
rz(-0.92241241) q[0];
x q[1];
rz(-1.0240779) q[2];
sx q[2];
rz(-0.77968979) q[2];
sx q[2];
rz(-2.0852154) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.071242407) q[1];
sx q[1];
rz(-1.6082941) q[1];
sx q[1];
rz(2.9686023) q[1];
x q[2];
rz(-2.6218518) q[3];
sx q[3];
rz(-1.9155353) q[3];
sx q[3];
rz(-2.910579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7543588) q[2];
sx q[2];
rz(-1.0789472) q[2];
sx q[2];
rz(-2.8299502) q[2];
rz(-1.257487) q[3];
sx q[3];
rz(-2.8897372) q[3];
sx q[3];
rz(-2.9529412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99609128) q[0];
sx q[0];
rz(-1.6956734) q[0];
sx q[0];
rz(-1.9022994) q[0];
rz(0.39341012) q[1];
sx q[1];
rz(-2.0937803) q[1];
sx q[1];
rz(-1.0308824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4811343) q[0];
sx q[0];
rz(-1.0054614) q[0];
sx q[0];
rz(-0.41445606) q[0];
x q[1];
rz(-1.8003045) q[2];
sx q[2];
rz(-1.3131427) q[2];
sx q[2];
rz(2.9837556) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46566612) q[1];
sx q[1];
rz(-0.74518004) q[1];
sx q[1];
rz(-2.096677) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83194701) q[3];
sx q[3];
rz(-2.0752677) q[3];
sx q[3];
rz(2.968806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3139412) q[2];
sx q[2];
rz(-0.88259077) q[2];
sx q[2];
rz(2.7871056) q[2];
rz(2.3101824) q[3];
sx q[3];
rz(-0.76787132) q[3];
sx q[3];
rz(-1.6262511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671234) q[0];
sx q[0];
rz(-2.9614083) q[0];
sx q[0];
rz(-2.6572976) q[0];
rz(-2.2593185) q[1];
sx q[1];
rz(-2.2703998) q[1];
sx q[1];
rz(-0.54214111) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6279468) q[0];
sx q[0];
rz(-1.8812669) q[0];
sx q[0];
rz(1.9438351) q[0];
rz(2.9464821) q[2];
sx q[2];
rz(-2.8144022) q[2];
sx q[2];
rz(-1.7649775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0391991) q[1];
sx q[1];
rz(-2.3619283) q[1];
sx q[1];
rz(-1.8236266) q[1];
rz(-1.3123973) q[3];
sx q[3];
rz(-0.2451788) q[3];
sx q[3];
rz(2.6441531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.813628) q[2];
sx q[2];
rz(-0.70983228) q[2];
sx q[2];
rz(2.1072809) q[2];
rz(1.3616925) q[3];
sx q[3];
rz(-1.9618278) q[3];
sx q[3];
rz(-2.5907607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66083241) q[0];
sx q[0];
rz(-1.3897422) q[0];
sx q[0];
rz(2.0939636) q[0];
rz(1.3230336) q[1];
sx q[1];
rz(-0.58224693) q[1];
sx q[1];
rz(0.38527647) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91222969) q[0];
sx q[0];
rz(-1.9989387) q[0];
sx q[0];
rz(1.0083126) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6462506) q[2];
sx q[2];
rz(-2.2860043) q[2];
sx q[2];
rz(-1.0714517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2400018) q[1];
sx q[1];
rz(-1.6419672) q[1];
sx q[1];
rz(-0.20889568) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10984595) q[3];
sx q[3];
rz(-0.89021909) q[3];
sx q[3];
rz(2.3504013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.72383991) q[2];
sx q[2];
rz(-2.1989792) q[2];
sx q[2];
rz(2.8670132) q[2];
rz(-2.5802021) q[3];
sx q[3];
rz(-1.0654819) q[3];
sx q[3];
rz(-1.5076465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1356337) q[0];
sx q[0];
rz(-0.57666403) q[0];
sx q[0];
rz(0.86719257) q[0];
rz(-2.7658956) q[1];
sx q[1];
rz(-0.58940327) q[1];
sx q[1];
rz(1.2295178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73342434) q[0];
sx q[0];
rz(-0.46981341) q[0];
sx q[0];
rz(0.3656268) q[0];
x q[1];
rz(-2.0784573) q[2];
sx q[2];
rz(-2.5637321) q[2];
sx q[2];
rz(-0.40921989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.83513488) q[1];
sx q[1];
rz(-2.0270438) q[1];
sx q[1];
rz(-0.92822509) q[1];
rz(-pi) q[2];
rz(0.8939871) q[3];
sx q[3];
rz(-0.5115307) q[3];
sx q[3];
rz(2.3681896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1162954) q[2];
sx q[2];
rz(-1.2325492) q[2];
sx q[2];
rz(-3.0787025) q[2];
rz(-2.7316015) q[3];
sx q[3];
rz(-2.4153109) q[3];
sx q[3];
rz(2.9819152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9922239) q[0];
sx q[0];
rz(-0.069644444) q[0];
sx q[0];
rz(0.83576354) q[0];
rz(1.958485) q[1];
sx q[1];
rz(-1.8712021) q[1];
sx q[1];
rz(-0.74367181) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27722699) q[0];
sx q[0];
rz(-2.1081308) q[0];
sx q[0];
rz(1.9242084) q[0];
rz(-pi) q[1];
rz(1.2116777) q[2];
sx q[2];
rz(-0.15378498) q[2];
sx q[2];
rz(-1.6395456) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3144296) q[1];
sx q[1];
rz(-1.3884228) q[1];
sx q[1];
rz(2.2807987) q[1];
rz(-pi) q[2];
rz(-1.623769) q[3];
sx q[3];
rz(-1.3154239) q[3];
sx q[3];
rz(2.6290174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30367294) q[2];
sx q[2];
rz(-0.49771365) q[2];
sx q[2];
rz(-0.79353235) q[2];
rz(0.41905904) q[3];
sx q[3];
rz(-1.4895118) q[3];
sx q[3];
rz(-0.52217531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94988743) q[0];
sx q[0];
rz(-0.97923034) q[0];
sx q[0];
rz(-1.1631843) q[0];
rz(-2.361182) q[1];
sx q[1];
rz(-2.8051832) q[1];
sx q[1];
rz(-1.6132678) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96161973) q[0];
sx q[0];
rz(-0.69937569) q[0];
sx q[0];
rz(-1.4699303) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0397908) q[2];
sx q[2];
rz(-2.7534979) q[2];
sx q[2];
rz(-2.3562252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0039499) q[1];
sx q[1];
rz(-0.19927916) q[1];
sx q[1];
rz(-1.6780705) q[1];
x q[2];
rz(2.1329857) q[3];
sx q[3];
rz(-2.6531124) q[3];
sx q[3];
rz(2.307596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10716001) q[2];
sx q[2];
rz(-1.9623423) q[2];
sx q[2];
rz(0.068664702) q[2];
rz(2.5717403) q[3];
sx q[3];
rz(-0.48044258) q[3];
sx q[3];
rz(-2.7710052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625075) q[0];
sx q[0];
rz(-2.469049) q[0];
sx q[0];
rz(1.8261209) q[0];
rz(-1.8585809) q[1];
sx q[1];
rz(-2.7048769) q[1];
sx q[1];
rz(0.049093094) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0822711) q[0];
sx q[0];
rz(-1.647864) q[0];
sx q[0];
rz(-1.661411) q[0];
x q[1];
rz(0.95048381) q[2];
sx q[2];
rz(-1.1178218) q[2];
sx q[2];
rz(0.17826232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0480501) q[1];
sx q[1];
rz(-0.88359264) q[1];
sx q[1];
rz(-0.50205135) q[1];
x q[2];
rz(-0.12709789) q[3];
sx q[3];
rz(-2.2107901) q[3];
sx q[3];
rz(2.618263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7580238) q[2];
sx q[2];
rz(-0.87791666) q[2];
sx q[2];
rz(2.6007268) q[2];
rz(1.0848378) q[3];
sx q[3];
rz(-2.4386051) q[3];
sx q[3];
rz(1.3802403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.908602) q[0];
sx q[0];
rz(-0.99642307) q[0];
sx q[0];
rz(-3.0365699) q[0];
rz(-2.5999293) q[1];
sx q[1];
rz(-0.88625208) q[1];
sx q[1];
rz(-2.7856316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5585238) q[0];
sx q[0];
rz(-1.5332869) q[0];
sx q[0];
rz(-1.7100699) q[0];
rz(-0.84635205) q[2];
sx q[2];
rz(-1.0484107) q[2];
sx q[2];
rz(-1.9047996) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4222504) q[1];
sx q[1];
rz(-1.0021035) q[1];
sx q[1];
rz(2.689792) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0186986) q[3];
sx q[3];
rz(-1.8884522) q[3];
sx q[3];
rz(1.7302707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2124704) q[2];
sx q[2];
rz(-1.6419623) q[2];
sx q[2];
rz(1.8222202) q[2];
rz(1.8170554) q[3];
sx q[3];
rz(-1.5732485) q[3];
sx q[3];
rz(0.96778473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20973715) q[0];
sx q[0];
rz(-3.1071438) q[0];
sx q[0];
rz(-1.6784278) q[0];
rz(1.6819008) q[1];
sx q[1];
rz(-1.1594783) q[1];
sx q[1];
rz(0.34585888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56014868) q[0];
sx q[0];
rz(-1.5236824) q[0];
sx q[0];
rz(-1.4909049) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5058636) q[2];
sx q[2];
rz(-1.254515) q[2];
sx q[2];
rz(-2.8313178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1766369) q[1];
sx q[1];
rz(-1.1208911) q[1];
sx q[1];
rz(0.69590203) q[1];
rz(-pi) q[2];
rz(2.8000051) q[3];
sx q[3];
rz(-2.1862967) q[3];
sx q[3];
rz(-0.08344354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2807002) q[2];
sx q[2];
rz(-1.8701376) q[2];
sx q[2];
rz(-0.26091179) q[2];
rz(-1.4704618) q[3];
sx q[3];
rz(-1.4025531) q[3];
sx q[3];
rz(-1.7117333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3093001) q[0];
sx q[0];
rz(-2.0335048) q[0];
sx q[0];
rz(-0.19620398) q[0];
rz(2.590754) q[1];
sx q[1];
rz(-1.486634) q[1];
sx q[1];
rz(-2.6015729) q[1];
rz(0.84340855) q[2];
sx q[2];
rz(-1.7512097) q[2];
sx q[2];
rz(-0.51633121) q[2];
rz(-2.4110386) q[3];
sx q[3];
rz(-1.3638221) q[3];
sx q[3];
rz(2.4331349) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
