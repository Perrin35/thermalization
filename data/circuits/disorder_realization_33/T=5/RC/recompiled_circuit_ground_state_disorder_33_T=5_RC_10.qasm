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
rz(9.8812048) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(2.1906617) q[1];
sx q[1];
rz(9.6074109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4066494) q[0];
sx q[0];
rz(-1.8514086) q[0];
sx q[0];
rz(-2.2191802) q[0];
rz(-pi) q[1];
rz(2.6668197) q[2];
sx q[2];
rz(-0.92657303) q[2];
sx q[2];
rz(-2.7930773) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4306972) q[1];
sx q[1];
rz(-0.17696807) q[1];
sx q[1];
rz(-2.9269993) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9630646) q[3];
sx q[3];
rz(-1.0844106) q[3];
sx q[3];
rz(1.148759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38723382) q[2];
sx q[2];
rz(-1.0789472) q[2];
sx q[2];
rz(-0.31164247) q[2];
rz(-1.8841057) q[3];
sx q[3];
rz(-2.8897372) q[3];
sx q[3];
rz(-0.18865147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99609128) q[0];
sx q[0];
rz(-1.6956734) q[0];
sx q[0];
rz(-1.2392932) q[0];
rz(-2.7481825) q[1];
sx q[1];
rz(-1.0478123) q[1];
sx q[1];
rz(-2.1107103) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14178947) q[0];
sx q[0];
rz(-1.9177525) q[0];
sx q[0];
rz(0.96472558) q[0];
x q[1];
rz(-2.4294063) q[2];
sx q[2];
rz(-2.7982494) q[2];
sx q[2];
rz(2.5573829) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0075025) q[1];
sx q[1];
rz(-2.1975027) q[1];
sx q[1];
rz(2.7078663) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64149789) q[3];
sx q[3];
rz(-0.94029762) q[3];
sx q[3];
rz(1.3288095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3139412) q[2];
sx q[2];
rz(-0.88259077) q[2];
sx q[2];
rz(2.7871056) q[2];
rz(-2.3101824) q[3];
sx q[3];
rz(-2.3737213) q[3];
sx q[3];
rz(-1.6262511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77446929) q[0];
sx q[0];
rz(-0.18018436) q[0];
sx q[0];
rz(-2.6572976) q[0];
rz(-0.88227415) q[1];
sx q[1];
rz(-0.87119281) q[1];
sx q[1];
rz(-0.54214111) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.079744) q[0];
sx q[0];
rz(-1.925191) q[0];
sx q[0];
rz(-0.3318048) q[0];
rz(0.32140478) q[2];
sx q[2];
rz(-1.5084477) q[2];
sx q[2];
rz(-0.37918175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24607813) q[1];
sx q[1];
rz(-2.3195004) q[1];
sx q[1];
rz(-2.8991634) q[1];
rz(-pi) q[2];
rz(-1.3334502) q[3];
sx q[3];
rz(-1.508731) q[3];
sx q[3];
rz(0.82236034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.813628) q[2];
sx q[2];
rz(-2.4317604) q[2];
sx q[2];
rz(1.0343118) q[2];
rz(-1.7799001) q[3];
sx q[3];
rz(-1.9618278) q[3];
sx q[3];
rz(-2.5907607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4807602) q[0];
sx q[0];
rz(-1.3897422) q[0];
sx q[0];
rz(2.0939636) q[0];
rz(1.3230336) q[1];
sx q[1];
rz(-0.58224693) q[1];
sx q[1];
rz(0.38527647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7390201) q[0];
sx q[0];
rz(-2.077335) q[0];
sx q[0];
rz(0.49474025) q[0];
rz(1.6462506) q[2];
sx q[2];
rz(-2.2860043) q[2];
sx q[2];
rz(-1.0714517) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99285728) q[1];
sx q[1];
rz(-0.22051935) q[1];
sx q[1];
rz(2.810477) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2543318) q[3];
sx q[3];
rz(-1.6561015) q[3];
sx q[3];
rz(2.4312756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72383991) q[2];
sx q[2];
rz(-2.1989792) q[2];
sx q[2];
rz(0.27457944) q[2];
rz(-0.56139055) q[3];
sx q[3];
rz(-1.0654819) q[3];
sx q[3];
rz(-1.6339462) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1356337) q[0];
sx q[0];
rz(-0.57666403) q[0];
sx q[0];
rz(-2.2744001) q[0];
rz(-0.37569702) q[1];
sx q[1];
rz(-2.5521894) q[1];
sx q[1];
rz(1.2295178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.00261) q[0];
sx q[0];
rz(-2.0073237) q[0];
sx q[0];
rz(-1.3912203) q[0];
rz(-pi) q[1];
rz(-2.0784573) q[2];
sx q[2];
rz(-0.57786059) q[2];
sx q[2];
rz(0.40921989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3064578) q[1];
sx q[1];
rz(-2.0270438) q[1];
sx q[1];
rz(-0.92822509) q[1];
rz(0.8939871) q[3];
sx q[3];
rz(-0.5115307) q[3];
sx q[3];
rz(2.3681896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1162954) q[2];
sx q[2];
rz(-1.9090434) q[2];
sx q[2];
rz(0.06289014) q[2];
rz(0.40999117) q[3];
sx q[3];
rz(-0.72628179) q[3];
sx q[3];
rz(-2.9819152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9922239) q[0];
sx q[0];
rz(-0.069644444) q[0];
sx q[0];
rz(-0.83576354) q[0];
rz(1.1831076) q[1];
sx q[1];
rz(-1.2703905) q[1];
sx q[1];
rz(-0.74367181) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8643657) q[0];
sx q[0];
rz(-1.0334618) q[0];
sx q[0];
rz(-1.2173843) q[0];
x q[1];
rz(-1.9299149) q[2];
sx q[2];
rz(-0.15378498) q[2];
sx q[2];
rz(-1.6395456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89827713) q[1];
sx q[1];
rz(-2.26663) q[1];
sx q[1];
rz(-0.23855539) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.623769) q[3];
sx q[3];
rz(-1.3154239) q[3];
sx q[3];
rz(-0.51257521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30367294) q[2];
sx q[2];
rz(-0.49771365) q[2];
sx q[2];
rz(0.79353235) q[2];
rz(-2.7225336) q[3];
sx q[3];
rz(-1.6520809) q[3];
sx q[3];
rz(-2.6194173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1917052) q[0];
sx q[0];
rz(-2.1623623) q[0];
sx q[0];
rz(1.9784084) q[0];
rz(2.361182) q[1];
sx q[1];
rz(-0.33640948) q[1];
sx q[1];
rz(-1.6132678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0484783) q[0];
sx q[0];
rz(-0.87568863) q[0];
sx q[0];
rz(3.0570875) q[0];
rz(-pi) q[1];
rz(-1.2210991) q[2];
sx q[2];
rz(-1.3989067) q[2];
sx q[2];
rz(0.34687172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1376428) q[1];
sx q[1];
rz(-0.19927916) q[1];
sx q[1];
rz(-1.6780705) q[1];
rz(-pi) q[2];
rz(-2.1329857) q[3];
sx q[3];
rz(-0.48848029) q[3];
sx q[3];
rz(2.307596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10716001) q[2];
sx q[2];
rz(-1.9623423) q[2];
sx q[2];
rz(-3.072928) q[2];
rz(-0.56985235) q[3];
sx q[3];
rz(-0.48044258) q[3];
sx q[3];
rz(0.3705875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625075) q[0];
sx q[0];
rz(-2.469049) q[0];
sx q[0];
rz(1.8261209) q[0];
rz(1.8585809) q[1];
sx q[1];
rz(-2.7048769) q[1];
sx q[1];
rz(-0.049093094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4955208) q[0];
sx q[0];
rz(-1.4804513) q[0];
sx q[0];
rz(-3.0642088) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2678864) q[2];
sx q[2];
rz(-0.75007909) q[2];
sx q[2];
rz(1.1994135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0480501) q[1];
sx q[1];
rz(-2.258) q[1];
sx q[1];
rz(2.6395413) q[1];
rz(-2.2146739) q[3];
sx q[3];
rz(-1.4689494) q[3];
sx q[3];
rz(-2.1702914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38356885) q[2];
sx q[2];
rz(-2.263676) q[2];
sx q[2];
rz(2.6007268) q[2];
rz(-2.0567549) q[3];
sx q[3];
rz(-2.4386051) q[3];
sx q[3];
rz(-1.7613523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.2553406) q[1];
sx q[1];
rz(-0.35596102) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5830688) q[0];
sx q[0];
rz(-1.5332869) q[0];
sx q[0];
rz(1.7100699) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2860724) q[2];
sx q[2];
rz(-0.86453712) q[2];
sx q[2];
rz(0.17937961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4527725) q[1];
sx q[1];
rz(-0.710462) q[1];
sx q[1];
rz(-0.97150357) q[1];
x q[2];
rz(-1.2139236) q[3];
sx q[3];
rz(-0.3398474) q[3];
sx q[3];
rz(1.7879037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2124704) q[2];
sx q[2];
rz(-1.6419623) q[2];
sx q[2];
rz(-1.8222202) q[2];
rz(-1.8170554) q[3];
sx q[3];
rz(-1.5732485) q[3];
sx q[3];
rz(-0.96778473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9318555) q[0];
sx q[0];
rz(-0.034448817) q[0];
sx q[0];
rz(1.4631648) q[0];
rz(1.4596918) q[1];
sx q[1];
rz(-1.1594783) q[1];
sx q[1];
rz(-0.34585888) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5992085) q[0];
sx q[0];
rz(-0.092723474) q[0];
sx q[0];
rz(-2.1044162) q[0];
rz(-pi) q[1];
rz(0.31690545) q[2];
sx q[2];
rz(-1.6325041) q[2];
sx q[2];
rz(1.8608492) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74259669) q[1];
sx q[1];
rz(-2.1861939) q[1];
sx q[1];
rz(-1.0091429) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92697619) q[3];
sx q[3];
rz(-1.8478353) q[3];
sx q[3];
rz(1.8566956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86089245) q[2];
sx q[2];
rz(-1.2714551) q[2];
sx q[2];
rz(-2.8806809) q[2];
rz(1.4704618) q[3];
sx q[3];
rz(-1.4025531) q[3];
sx q[3];
rz(-1.4298593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83229257) q[0];
sx q[0];
rz(-1.1080879) q[0];
sx q[0];
rz(2.9453887) q[0];
rz(2.590754) q[1];
sx q[1];
rz(-1.486634) q[1];
sx q[1];
rz(-2.6015729) q[1];
rz(0.2395128) q[2];
sx q[2];
rz(-2.2838267) q[2];
sx q[2];
rz(-2.2455278) q[2];
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
