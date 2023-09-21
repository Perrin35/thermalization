OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(4.5933525) q[0];
sx q[0];
rz(10.070355) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49255532) q[0];
sx q[0];
rz(-1.4264002) q[0];
sx q[0];
rz(1.210968) q[0];
rz(-pi) q[1];
rz(3.0718578) q[2];
sx q[2];
rz(-2.2001534) q[2];
sx q[2];
rz(-2.8024408) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9420535) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(-2.8295423) q[1];
rz(-pi) q[2];
rz(1.9707768) q[3];
sx q[3];
rz(-2.1809289) q[3];
sx q[3];
rz(-0.67584544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-2.5966068) q[0];
rz(0.90822059) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(-0.82495904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075515) q[0];
sx q[0];
rz(-1.3836765) q[0];
sx q[0];
rz(1.1603111) q[0];
rz(3.1171563) q[2];
sx q[2];
rz(-2.7385798) q[2];
sx q[2];
rz(-1.8624061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7823239) q[1];
sx q[1];
rz(-1.4546848) q[1];
sx q[1];
rz(-2.5776723) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0042079) q[3];
sx q[3];
rz(-1.7363747) q[3];
sx q[3];
rz(-2.2752938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58723441) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(-0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3787057) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(3.1030531) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24667106) q[0];
sx q[0];
rz(-1.5630432) q[0];
sx q[0];
rz(0.00064108032) q[0];
rz(-pi) q[1];
rz(-1.4171717) q[2];
sx q[2];
rz(-1.7989752) q[2];
sx q[2];
rz(1.7287901) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.027027834) q[1];
sx q[1];
rz(-1.6323043) q[1];
sx q[1];
rz(-2.0608749) q[1];
x q[2];
rz(0.97088082) q[3];
sx q[3];
rz(-2.0754793) q[3];
sx q[3];
rz(1.8463299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58275756) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38347605) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-0.19515881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.450481) q[0];
sx q[0];
rz(-3.0483732) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0267369) q[2];
sx q[2];
rz(-1.1508905) q[2];
sx q[2];
rz(-1.3865711) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2745167) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(2.0348674) q[1];
rz(0.021899453) q[3];
sx q[3];
rz(-2.5023513) q[3];
sx q[3];
rz(-2.1454449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347663) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(-1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(2.7635014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3867214) q[0];
sx q[0];
rz(-0.62749642) q[0];
sx q[0];
rz(2.0621215) q[0];
rz(2.3230882) q[2];
sx q[2];
rz(-0.92220014) q[2];
sx q[2];
rz(-1.6550145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8564344) q[1];
sx q[1];
rz(-2.563623) q[1];
sx q[1];
rz(-2.2846089) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37851815) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(-0.62687031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(0.88796973) q[2];
rz(-0.97638431) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-0.37762541) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(0.91517085) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5765502) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(-2.8898426) q[0];
rz(-pi) q[1];
rz(2.633811) q[2];
sx q[2];
rz(-1.0630597) q[2];
sx q[2];
rz(2.7024262) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8296216) q[1];
sx q[1];
rz(-1.8835856) q[1];
sx q[1];
rz(-1.4986067) q[1];
rz(0.69453199) q[3];
sx q[3];
rz(-1.4834705) q[3];
sx q[3];
rz(-0.15184034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(0.79052314) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054366) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(-0.36703584) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050497748) q[0];
sx q[0];
rz(-2.0031345) q[0];
sx q[0];
rz(1.7763441) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13397127) q[2];
sx q[2];
rz(-0.44951648) q[2];
sx q[2];
rz(-1.9884895) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62896171) q[1];
sx q[1];
rz(-0.96620622) q[1];
sx q[1];
rz(1.1058034) q[1];
rz(-pi) q[2];
rz(-0.013399259) q[3];
sx q[3];
rz(-1.0671167) q[3];
sx q[3];
rz(-0.79813938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(-1.210775) q[2];
rz(0.0018421729) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84412557) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-2.3892367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90826666) q[0];
sx q[0];
rz(-1.255863) q[0];
sx q[0];
rz(0.53091913) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3812177) q[2];
sx q[2];
rz(-0.65738064) q[2];
sx q[2];
rz(-1.1084523) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89531089) q[1];
sx q[1];
rz(-1.6591676) q[1];
sx q[1];
rz(-1.4573216) q[1];
rz(-0.30267834) q[3];
sx q[3];
rz(-0.61121002) q[3];
sx q[3];
rz(-1.4509033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2032808) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-3.1414462) q[2];
rz(1.1095307) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-1.0060271) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.8267652) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6515503) q[0];
sx q[0];
rz(-2.393912) q[0];
sx q[0];
rz(2.7689395) q[0];
x q[1];
rz(-1.7039653) q[2];
sx q[2];
rz(-1.9317992) q[2];
sx q[2];
rz(-1.5683057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4970376) q[1];
sx q[1];
rz(-2.2215448) q[1];
sx q[1];
rz(0.52849309) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5789088) q[3];
sx q[3];
rz(-2.3981961) q[3];
sx q[3];
rz(-2.5826366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7018147) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(-0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4196639) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(-0.28668177) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-0.39189664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4927917) q[0];
sx q[0];
rz(-1.0133044) q[0];
sx q[0];
rz(-3.1317943) q[0];
x q[1];
rz(1.7461807) q[2];
sx q[2];
rz(-0.66329623) q[2];
sx q[2];
rz(-1.7593918) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0185768) q[1];
sx q[1];
rz(-1.3850817) q[1];
sx q[1];
rz(2.3260444) q[1];
rz(-pi) q[2];
rz(2.6565353) q[3];
sx q[3];
rz(-2.5254446) q[3];
sx q[3];
rz(-0.96998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60951704) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(1.9840476) q[2];
rz(2.5907497) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75523238) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(1.8023087) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(2.9677283) q[2];
sx q[2];
rz(-2.3542913) q[2];
sx q[2];
rz(-2.3402294) q[2];
rz(-0.34655456) q[3];
sx q[3];
rz(-2.0934436) q[3];
sx q[3];
rz(2.2128076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
