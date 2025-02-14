OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27788568) q[0];
sx q[0];
rz(-2.7122893) q[0];
sx q[0];
rz(2.8306146) q[0];
rz(-1.3485981) q[1];
sx q[1];
rz(4.2378874) q[1];
sx q[1];
rz(8.8506946) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1659083) q[0];
sx q[0];
rz(-0.57281369) q[0];
sx q[0];
rz(2.2982135) q[0];
rz(-pi) q[1];
rz(1.6234196) q[2];
sx q[2];
rz(-0.5835909) q[2];
sx q[2];
rz(2.6578195) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8694386) q[1];
sx q[1];
rz(-1.660214) q[1];
sx q[1];
rz(-1.2050259) q[1];
x q[2];
rz(-1.7314243) q[3];
sx q[3];
rz(-2.7607252) q[3];
sx q[3];
rz(-0.67202744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2744039) q[2];
sx q[2];
rz(-1.8827266) q[2];
sx q[2];
rz(1.8223507) q[2];
rz(0.073171767) q[3];
sx q[3];
rz(-1.3061482) q[3];
sx q[3];
rz(-0.81899548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0935593) q[0];
sx q[0];
rz(-1.9608542) q[0];
sx q[0];
rz(-2.4617885) q[0];
rz(-2.3381084) q[1];
sx q[1];
rz(-1.7796703) q[1];
sx q[1];
rz(-0.63978535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6006192) q[0];
sx q[0];
rz(-2.2007211) q[0];
sx q[0];
rz(0.39430228) q[0];
rz(-1.4714986) q[2];
sx q[2];
rz(-2.298215) q[2];
sx q[2];
rz(-2.119144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27187335) q[1];
sx q[1];
rz(-0.98733989) q[1];
sx q[1];
rz(2.3223552) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62044797) q[3];
sx q[3];
rz(-0.37319601) q[3];
sx q[3];
rz(2.9377112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.60376072) q[2];
sx q[2];
rz(-0.49850264) q[2];
sx q[2];
rz(-2.4375088) q[2];
rz(-0.1121366) q[3];
sx q[3];
rz(-1.6928558) q[3];
sx q[3];
rz(1.0224379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10244399) q[0];
sx q[0];
rz(-1.1071438) q[0];
sx q[0];
rz(0.33945864) q[0];
rz(2.0231694) q[1];
sx q[1];
rz(-2.2148841) q[1];
sx q[1];
rz(3.1059473) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46860002) q[0];
sx q[0];
rz(-2.9170999) q[0];
sx q[0];
rz(-0.39546449) q[0];
rz(1.8730803) q[2];
sx q[2];
rz(-1.8955823) q[2];
sx q[2];
rz(-3.0168282) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0633723) q[1];
sx q[1];
rz(-0.92670346) q[1];
sx q[1];
rz(-1.3011342) q[1];
rz(0.6142637) q[3];
sx q[3];
rz(-1.302056) q[3];
sx q[3];
rz(-0.72878557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3493335) q[2];
sx q[2];
rz(-2.4534241) q[2];
sx q[2];
rz(-2.2815857) q[2];
rz(-0.79259121) q[3];
sx q[3];
rz(-0.20321295) q[3];
sx q[3];
rz(-2.121076) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0122796) q[0];
sx q[0];
rz(-1.3849994) q[0];
sx q[0];
rz(-2.5087575) q[0];
rz(-1.1526147) q[1];
sx q[1];
rz(-0.8747789) q[1];
sx q[1];
rz(-1.4702183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9235719) q[0];
sx q[0];
rz(-1.583719) q[0];
sx q[0];
rz(0.92148975) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7196413) q[2];
sx q[2];
rz(-0.80721569) q[2];
sx q[2];
rz(2.5850353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51293514) q[1];
sx q[1];
rz(-2.3875065) q[1];
sx q[1];
rz(-1.1910466) q[1];
x q[2];
rz(2.4864462) q[3];
sx q[3];
rz(-2.0787424) q[3];
sx q[3];
rz(0.20482132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4319438) q[2];
sx q[2];
rz(-2.180763) q[2];
sx q[2];
rz(-2.7107837) q[2];
rz(0.68239799) q[3];
sx q[3];
rz(-1.4902481) q[3];
sx q[3];
rz(0.94990134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7948941) q[0];
sx q[0];
rz(-1.8759202) q[0];
sx q[0];
rz(-1.1899813) q[0];
rz(-0.31040141) q[1];
sx q[1];
rz(-2.0605395) q[1];
sx q[1];
rz(0.50500542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0520262) q[0];
sx q[0];
rz(-2.1329693) q[0];
sx q[0];
rz(-0.52345353) q[0];
rz(-pi) q[1];
rz(1.861635) q[2];
sx q[2];
rz(-2.7105936) q[2];
sx q[2];
rz(0.19311007) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6663899) q[1];
sx q[1];
rz(-1.0745) q[1];
sx q[1];
rz(-1.6700891) q[1];
x q[2];
rz(-0.72568958) q[3];
sx q[3];
rz(-1.6349051) q[3];
sx q[3];
rz(0.025246092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7348822) q[2];
sx q[2];
rz(-0.70047417) q[2];
sx q[2];
rz(0.90275466) q[2];
rz(2.086575) q[3];
sx q[3];
rz(-2.2746634) q[3];
sx q[3];
rz(-2.6796851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2257253) q[0];
sx q[0];
rz(-2.4287455) q[0];
sx q[0];
rz(-1.073904) q[0];
rz(-1.3794927) q[1];
sx q[1];
rz(-2.4122489) q[1];
sx q[1];
rz(3.0440547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65501761) q[0];
sx q[0];
rz(-1.7518105) q[0];
sx q[0];
rz(0.54365309) q[0];
x q[1];
rz(-1.9768841) q[2];
sx q[2];
rz(-2.7345719) q[2];
sx q[2];
rz(-0.87813745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7382724) q[1];
sx q[1];
rz(-3.0636859) q[1];
sx q[1];
rz(2.2408443) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0076526) q[3];
sx q[3];
rz(-1.1112257) q[3];
sx q[3];
rz(-0.28467049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11233106) q[2];
sx q[2];
rz(-1.4533726) q[2];
sx q[2];
rz(0.41109273) q[2];
rz(-1.4136275) q[3];
sx q[3];
rz(-2.3634383) q[3];
sx q[3];
rz(-1.5948064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77808648) q[0];
sx q[0];
rz(-0.030487617) q[0];
sx q[0];
rz(1.4084858) q[0];
rz(-2.1259437) q[1];
sx q[1];
rz(-1.8135095) q[1];
sx q[1];
rz(-1.4782864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41970872) q[0];
sx q[0];
rz(-0.96707487) q[0];
sx q[0];
rz(1.2092071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99009606) q[2];
sx q[2];
rz(-2.1352571) q[2];
sx q[2];
rz(1.1549468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6802753) q[1];
sx q[1];
rz(-1.6467386) q[1];
sx q[1];
rz(-1.126775) q[1];
x q[2];
rz(1.8084548) q[3];
sx q[3];
rz(-2.6293652) q[3];
sx q[3];
rz(0.32853261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0279072) q[2];
sx q[2];
rz(-2.5137641) q[2];
sx q[2];
rz(0.17024635) q[2];
rz(-1.8611192) q[3];
sx q[3];
rz(-1.6672641) q[3];
sx q[3];
rz(1.1580275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.668642) q[0];
sx q[0];
rz(-0.96422115) q[0];
sx q[0];
rz(-0.52841312) q[0];
rz(-0.93327418) q[1];
sx q[1];
rz(-0.92527881) q[1];
sx q[1];
rz(-2.5859213) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816102) q[0];
sx q[0];
rz(-0.9564119) q[0];
sx q[0];
rz(-0.23991628) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.470475) q[2];
sx q[2];
rz(-2.360376) q[2];
sx q[2];
rz(0.90512102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56951504) q[1];
sx q[1];
rz(-1.6028499) q[1];
sx q[1];
rz(-0.7116913) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8818086) q[3];
sx q[3];
rz(-2.2385983) q[3];
sx q[3];
rz(3.1329114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8133424) q[2];
sx q[2];
rz(-0.80208653) q[2];
sx q[2];
rz(0.30725202) q[2];
rz(-0.79636374) q[3];
sx q[3];
rz(-2.4107404) q[3];
sx q[3];
rz(-0.097631924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55427134) q[0];
sx q[0];
rz(-1.048943) q[0];
sx q[0];
rz(1.3294504) q[0];
rz(-0.21996552) q[1];
sx q[1];
rz(-0.34383067) q[1];
sx q[1];
rz(1.8230009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23585478) q[0];
sx q[0];
rz(-1.3038583) q[0];
sx q[0];
rz(-0.22375317) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4542142) q[2];
sx q[2];
rz(-1.0739084) q[2];
sx q[2];
rz(-1.0036482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4798051) q[1];
sx q[1];
rz(-1.1068245) q[1];
sx q[1];
rz(-2.7590092) q[1];
rz(-pi) q[2];
rz(-2.044146) q[3];
sx q[3];
rz(-2.2339377) q[3];
sx q[3];
rz(1.9846265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69507504) q[2];
sx q[2];
rz(-1.7138285) q[2];
sx q[2];
rz(1.2850777) q[2];
rz(2.255693) q[3];
sx q[3];
rz(-1.3934087) q[3];
sx q[3];
rz(1.5281965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1000243) q[0];
sx q[0];
rz(-0.73606857) q[0];
sx q[0];
rz(-2.8247483) q[0];
rz(0.44081229) q[1];
sx q[1];
rz(-2.3471954) q[1];
sx q[1];
rz(-2.7712834) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0531385) q[0];
sx q[0];
rz(-1.1897108) q[0];
sx q[0];
rz(0.033381406) q[0];
x q[1];
rz(-0.65014) q[2];
sx q[2];
rz(-0.95528379) q[2];
sx q[2];
rz(0.56734771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1052484) q[1];
sx q[1];
rz(-0.40648983) q[1];
sx q[1];
rz(2.2986733) q[1];
rz(-0.050906128) q[3];
sx q[3];
rz(-1.6744782) q[3];
sx q[3];
rz(-3.0391676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73654282) q[2];
sx q[2];
rz(-2.6506347) q[2];
sx q[2];
rz(1.5092108) q[2];
rz(2.3368808) q[3];
sx q[3];
rz(-1.5038306) q[3];
sx q[3];
rz(0.10778431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0263154) q[0];
sx q[0];
rz(-1.4125217) q[0];
sx q[0];
rz(-1.9052292) q[0];
rz(-0.23030494) q[1];
sx q[1];
rz(-2.8254012) q[1];
sx q[1];
rz(0.91513035) q[1];
rz(0.66441734) q[2];
sx q[2];
rz(-1.0894486) q[2];
sx q[2];
rz(1.1846381) q[2];
rz(0.32634278) q[3];
sx q[3];
rz(-0.97728609) q[3];
sx q[3];
rz(1.4090007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
