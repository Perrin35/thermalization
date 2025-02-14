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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(-0.43098488) q[0];
rz(0.35297901) q[1];
sx q[1];
rz(-1.9184435) q[1];
sx q[1];
rz(-2.7117742) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0551222) q[0];
sx q[0];
rz(-1.6608547) q[0];
sx q[0];
rz(-0.49205972) q[0];
rz(-pi) q[1];
rz(-2.0031646) q[2];
sx q[2];
rz(-1.585344) q[2];
sx q[2];
rz(-2.5231139) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1466521) q[1];
sx q[1];
rz(-1.2122248) q[1];
sx q[1];
rz(1.7943804) q[1];
x q[2];
rz(2.6453703) q[3];
sx q[3];
rz(-1.2321951) q[3];
sx q[3];
rz(-0.34211516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.063804403) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(1.0112313) q[2];
rz(-2.1003335) q[3];
sx q[3];
rz(-2.6285089) q[3];
sx q[3];
rz(0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0332727) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(-0.74747768) q[0];
rz(0.29391089) q[1];
sx q[1];
rz(-2.0103318) q[1];
sx q[1];
rz(0.25516137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19758148) q[0];
sx q[0];
rz(-1.4048049) q[0];
sx q[0];
rz(2.152488) q[0];
rz(-pi) q[1];
rz(-2.8361144) q[2];
sx q[2];
rz(-1.7490088) q[2];
sx q[2];
rz(-0.62159789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9772759) q[1];
sx q[1];
rz(-2.212388) q[1];
sx q[1];
rz(2.5733828) q[1];
x q[2];
rz(-0.94631715) q[3];
sx q[3];
rz(-1.4186275) q[3];
sx q[3];
rz(0.63380235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.030563844) q[2];
sx q[2];
rz(-0.054628987) q[2];
sx q[2];
rz(0.89152208) q[2];
rz(-0.24957481) q[3];
sx q[3];
rz(-2.2587743) q[3];
sx q[3];
rz(-2.7454624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7340649) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(-0.044128142) q[0];
rz(1.2491501) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(0.485802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7120948) q[0];
sx q[0];
rz(-2.9154239) q[0];
sx q[0];
rz(-0.087271376) q[0];
rz(-pi) q[1];
rz(-2.2052231) q[2];
sx q[2];
rz(-2.9108725) q[2];
sx q[2];
rz(-1.3619193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9618157) q[1];
sx q[1];
rz(-1.0423585) q[1];
sx q[1];
rz(-0.68317271) q[1];
x q[2];
rz(0.39704571) q[3];
sx q[3];
rz(-1.9730803) q[3];
sx q[3];
rz(-2.8232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8140063) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(-0.69474727) q[2];
rz(-0.57573777) q[3];
sx q[3];
rz(-1.271194) q[3];
sx q[3];
rz(-1.7339138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95076743) q[0];
sx q[0];
rz(-2.8468843) q[0];
sx q[0];
rz(-0.29275352) q[0];
rz(2.5726908) q[1];
sx q[1];
rz(-2.3141373) q[1];
sx q[1];
rz(-0.84691602) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50668651) q[0];
sx q[0];
rz(-1.4932079) q[0];
sx q[0];
rz(-1.1617178) q[0];
rz(-3.0693502) q[2];
sx q[2];
rz(-2.2390167) q[2];
sx q[2];
rz(0.39230686) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9255815) q[1];
sx q[1];
rz(-2.8329509) q[1];
sx q[1];
rz(0.68283783) q[1];
x q[2];
rz(-1.2388703) q[3];
sx q[3];
rz(-2.9369825) q[3];
sx q[3];
rz(-1.0910891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7920821) q[2];
sx q[2];
rz(-1.258472) q[2];
sx q[2];
rz(-2.9328031) q[2];
rz(-0.71792349) q[3];
sx q[3];
rz(-2.8585298) q[3];
sx q[3];
rz(0.94055241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013041) q[0];
sx q[0];
rz(-2.6191819) q[0];
sx q[0];
rz(2.8029602) q[0];
rz(-0.70308095) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(2.3032761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.188952) q[0];
sx q[0];
rz(-1.3172272) q[0];
sx q[0];
rz(2.6978561) q[0];
rz(-pi) q[1];
rz(0.19836004) q[2];
sx q[2];
rz(-0.83399978) q[2];
sx q[2];
rz(-1.5892346) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5456713) q[1];
sx q[1];
rz(-1.6041846) q[1];
sx q[1];
rz(2.2708793) q[1];
rz(-pi) q[2];
rz(0.57788604) q[3];
sx q[3];
rz(-2.4494723) q[3];
sx q[3];
rz(-0.96847615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4579953) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(2.1387157) q[2];
rz(-2.9966127) q[3];
sx q[3];
rz(-0.093791157) q[3];
sx q[3];
rz(0.063044757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1231287) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(0.96328324) q[0];
rz(-2.9604498) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(1.9021665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9299663) q[0];
sx q[0];
rz(-2.3986882) q[0];
sx q[0];
rz(2.5725098) q[0];
x q[1];
rz(0.507429) q[2];
sx q[2];
rz(-0.52971887) q[2];
sx q[2];
rz(1.5629753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9628004) q[1];
sx q[1];
rz(-1.9427248) q[1];
sx q[1];
rz(-0.40403251) q[1];
rz(1.419098) q[3];
sx q[3];
rz(-1.6983508) q[3];
sx q[3];
rz(-0.90326819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2565101) q[2];
sx q[2];
rz(-0.91826597) q[2];
sx q[2];
rz(0.93290848) q[2];
rz(1.4126623) q[3];
sx q[3];
rz(-1.7191073) q[3];
sx q[3];
rz(0.2093813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6949718) q[0];
sx q[0];
rz(-2.7734723) q[0];
sx q[0];
rz(2.3387261) q[0];
rz(-2.39957) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-2.7313357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78722505) q[0];
sx q[0];
rz(-2.4929308) q[0];
sx q[0];
rz(2.5543268) q[0];
rz(-pi) q[1];
x q[1];
rz(1.688399) q[2];
sx q[2];
rz(-0.80999331) q[2];
sx q[2];
rz(1.2076898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8150543) q[1];
sx q[1];
rz(-1.6449303) q[1];
sx q[1];
rz(1.733041) q[1];
rz(-2.5531261) q[3];
sx q[3];
rz(-2.3938673) q[3];
sx q[3];
rz(2.6315053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.86847574) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(-2.9203019) q[2];
rz(2.7224329) q[3];
sx q[3];
rz(-2.6594888) q[3];
sx q[3];
rz(-0.47590772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047423) q[0];
sx q[0];
rz(-0.99822799) q[0];
sx q[0];
rz(-1.5297484) q[0];
rz(1.3329685) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(-1.6882247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3842363) q[0];
sx q[0];
rz(-2.1290551) q[0];
sx q[0];
rz(0.7406969) q[0];
x q[1];
rz(1.6906934) q[2];
sx q[2];
rz(-1.3812307) q[2];
sx q[2];
rz(1.3055064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.61078) q[1];
sx q[1];
rz(-1.7208092) q[1];
sx q[1];
rz(-0.71041815) q[1];
rz(-pi) q[2];
rz(-1.8492886) q[3];
sx q[3];
rz(-1.4208111) q[3];
sx q[3];
rz(1.5453135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86614418) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(-0.010738372) q[2];
rz(-2.8209316) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(-0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7084932) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(-0.47472111) q[0];
rz(-2.8482598) q[1];
sx q[1];
rz(-2.0359998) q[1];
sx q[1];
rz(1.4795823) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2441263) q[0];
sx q[0];
rz(-3.0871824) q[0];
sx q[0];
rz(-2.3136086) q[0];
rz(-0.29724289) q[2];
sx q[2];
rz(-0.69801509) q[2];
sx q[2];
rz(0.4618123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.41614) q[1];
sx q[1];
rz(-1.801277) q[1];
sx q[1];
rz(-0.47941717) q[1];
rz(-pi) q[2];
rz(-0.94274272) q[3];
sx q[3];
rz(-2.8398879) q[3];
sx q[3];
rz(2.7127271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.097229615) q[2];
sx q[2];
rz(-1.2791415) q[2];
sx q[2];
rz(-0.43107671) q[2];
rz(-0.19032446) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(-2.2569807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97656074) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(-0.2188368) q[0];
rz(-0.95895514) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(-1.7494019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9746124) q[0];
sx q[0];
rz(-1.3777527) q[0];
sx q[0];
rz(1.6231322) q[0];
x q[1];
rz(-1.7151681) q[2];
sx q[2];
rz(-1.1544268) q[2];
sx q[2];
rz(1.8954111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47962418) q[1];
sx q[1];
rz(-1.2733409) q[1];
sx q[1];
rz(0.65315078) q[1];
rz(2.1851319) q[3];
sx q[3];
rz(-1.9298565) q[3];
sx q[3];
rz(-0.39310716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2739111) q[2];
sx q[2];
rz(-0.76649222) q[2];
sx q[2];
rz(-1.319818) q[2];
rz(-2.6093318) q[3];
sx q[3];
rz(-1.3479193) q[3];
sx q[3];
rz(1.5925315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821447) q[0];
sx q[0];
rz(-1.5503217) q[0];
sx q[0];
rz(0.4009608) q[0];
rz(0.22790146) q[1];
sx q[1];
rz(-2.0961998) q[1];
sx q[1];
rz(-1.6963522) q[1];
rz(2.5980742) q[2];
sx q[2];
rz(-2.666011) q[2];
sx q[2];
rz(0.70071349) q[2];
rz(0.99146546) q[3];
sx q[3];
rz(-1.4351063) q[3];
sx q[3];
rz(-1.2096005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
