OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(-0.0052069081) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(4.2978573) q[1];
sx q[1];
rz(8.2351091) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.386728) q[0];
sx q[0];
rz(-1.0611738) q[0];
sx q[0];
rz(1.9652363) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18007937) q[2];
sx q[2];
rz(-1.4305978) q[2];
sx q[2];
rz(0.3917429) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7835044) q[1];
sx q[1];
rz(-1.7184098) q[1];
sx q[1];
rz(1.258177) q[1];
rz(-pi) q[2];
rz(0.94250836) q[3];
sx q[3];
rz(-1.2357467) q[3];
sx q[3];
rz(3.0818919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(2.544196) q[2];
rz(-1.3655837) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269203) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(3.0283668) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118606) q[0];
sx q[0];
rz(-1.1456053) q[0];
sx q[0];
rz(-0.042516275) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8791061) q[2];
sx q[2];
rz(-3.065735) q[2];
sx q[2];
rz(-2.5469317) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7659521) q[1];
sx q[1];
rz(-1.1728219) q[1];
sx q[1];
rz(2.6938733) q[1];
rz(-pi) q[2];
x q[2];
rz(1.232997) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(-0.88463569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-2.9193027) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(-2.0942028) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(-3.0139794) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022284431) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(-0.22247252) q[0];
x q[1];
rz(-2.3348546) q[2];
sx q[2];
rz(-2.2819937) q[2];
sx q[2];
rz(-2.9850609) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4096654) q[1];
sx q[1];
rz(-1.8095784) q[1];
sx q[1];
rz(3.0569397) q[1];
rz(-pi) q[2];
rz(-2.8887799) q[3];
sx q[3];
rz(-2.1826934) q[3];
sx q[3];
rz(1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7814653) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(-2.8313417) q[2];
rz(-0.34960738) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629102) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(2.7754916) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(-2.7691832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.303064) q[0];
sx q[0];
rz(-1.6922564) q[0];
sx q[0];
rz(-2.2785506) q[0];
x q[1];
rz(0.98948688) q[2];
sx q[2];
rz(-1.2876236) q[2];
sx q[2];
rz(-2.1894933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6652674) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(-1.9515576) q[1];
rz(-pi) q[2];
rz(0.59471547) q[3];
sx q[3];
rz(-2.3200431) q[3];
sx q[3];
rz(-2.6298414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7238414) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(-0.27238971) q[2];
rz(0.90879905) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(-2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(-2.0102665) q[0];
rz(2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(2.4647443) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9899321) q[0];
sx q[0];
rz(-0.44758546) q[0];
sx q[0];
rz(0.18330343) q[0];
x q[1];
rz(2.4104775) q[2];
sx q[2];
rz(-1.7321246) q[2];
sx q[2];
rz(-2.4246755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10830282) q[1];
sx q[1];
rz(-1.8255594) q[1];
sx q[1];
rz(1.5441896) q[1];
rz(-0.058651794) q[3];
sx q[3];
rz(-0.39415112) q[3];
sx q[3];
rz(-2.4413721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0040434917) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.1452902) q[2];
rz(-0.14906135) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(1.0173652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260219) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(-2.3964264) q[1];
sx q[1];
rz(-2.7909653) q[1];
sx q[1];
rz(0.93313342) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60676735) q[0];
sx q[0];
rz(-1.2762428) q[0];
sx q[0];
rz(2.0539001) q[0];
x q[1];
rz(-2.6283693) q[2];
sx q[2];
rz(-0.66868082) q[2];
sx q[2];
rz(0.19323397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3999403) q[1];
sx q[1];
rz(-1.5730904) q[1];
sx q[1];
rz(-1.9386577) q[1];
rz(-0.17721456) q[3];
sx q[3];
rz(-0.92068499) q[3];
sx q[3];
rz(-0.23055102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(-0.28665001) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(-1.6749143) q[0];
rz(-2.1215227) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(1.0010304) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9065549) q[0];
sx q[0];
rz(-2.8585003) q[0];
sx q[0];
rz(-1.4366158) q[0];
rz(-pi) q[1];
rz(0.34044388) q[2];
sx q[2];
rz(-1.2203487) q[2];
sx q[2];
rz(0.0705138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.475358) q[1];
sx q[1];
rz(-1.959414) q[1];
sx q[1];
rz(2.9734441) q[1];
rz(-pi) q[2];
rz(2.4909199) q[3];
sx q[3];
rz(-1.7400029) q[3];
sx q[3];
rz(2.8229439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(0.10350791) q[2];
rz(-2.5497656) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(-0.6119588) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-2.7517095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.10605) q[0];
sx q[0];
rz(-0.70362008) q[0];
sx q[0];
rz(-2.5662867) q[0];
x q[1];
rz(0.50750081) q[2];
sx q[2];
rz(-1.9365053) q[2];
sx q[2];
rz(-2.4909004) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.08323616) q[1];
sx q[1];
rz(-0.55587686) q[1];
sx q[1];
rz(-2.2298989) q[1];
x q[2];
rz(-0.26384683) q[3];
sx q[3];
rz(-0.90875328) q[3];
sx q[3];
rz(2.2784233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-0.71425444) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(2.3642448) q[0];
rz(-2.2946987) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-0.27639595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696635) q[0];
sx q[0];
rz(-1.164325) q[0];
sx q[0];
rz(1.7448234) q[0];
rz(1.3771636) q[2];
sx q[2];
rz(-1.3467222) q[2];
sx q[2];
rz(-2.5335238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.285187) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(1.3473947) q[1];
rz(-pi) q[2];
rz(1.3235839) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(-2.0127399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(-3.1048807) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(-1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(-2.1113077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093826483) q[0];
sx q[0];
rz(-2.2802135) q[0];
sx q[0];
rz(-1.312027) q[0];
rz(-pi) q[1];
rz(1.2376386) q[2];
sx q[2];
rz(-0.91234708) q[2];
sx q[2];
rz(-1.2308987) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2620169) q[1];
sx q[1];
rz(-2.587662) q[1];
sx q[1];
rz(-2.5074156) q[1];
rz(-pi) q[2];
rz(2.685931) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(1.1841707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-0.34118787) q[2];
rz(-2.4009005) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(-0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13070233) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(2.4717992) q[1];
sx q[1];
rz(-2.5302946) q[1];
sx q[1];
rz(-1.5763462) q[1];
rz(1.3763935) q[2];
sx q[2];
rz(-1.1444848) q[2];
sx q[2];
rz(-0.3111006) q[2];
rz(-0.68941152) q[3];
sx q[3];
rz(-1.1082311) q[3];
sx q[3];
rz(-2.5267596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
