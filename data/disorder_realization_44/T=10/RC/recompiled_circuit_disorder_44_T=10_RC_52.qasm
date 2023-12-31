OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(3.0043998) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(2.611673) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7007028) q[0];
sx q[0];
rz(-0.39632495) q[0];
sx q[0];
rz(-1.8519782) q[0];
rz(0.86962236) q[2];
sx q[2];
rz(-2.6343971) q[2];
sx q[2];
rz(1.5324355) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9143608) q[1];
sx q[1];
rz(-2.025361) q[1];
sx q[1];
rz(-2.582002) q[1];
rz(-pi) q[2];
rz(-2.8785273) q[3];
sx q[3];
rz(-0.77583003) q[3];
sx q[3];
rz(-2.8924931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(2.8049862) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(-1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34823725) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(1.1938098) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(2.3056727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7556691) q[0];
sx q[0];
rz(-0.14005157) q[0];
sx q[0];
rz(-1.7007909) q[0];
x q[1];
rz(0.95401986) q[2];
sx q[2];
rz(-2.4980133) q[2];
sx q[2];
rz(-1.8804903) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.086120124) q[1];
sx q[1];
rz(-1.2985843) q[1];
sx q[1];
rz(-2.230456) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73987506) q[3];
sx q[3];
rz(-1.2736819) q[3];
sx q[3];
rz(1.5996931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8643643) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(-0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42831746) q[0];
sx q[0];
rz(-1.0773766) q[0];
sx q[0];
rz(2.0879478) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-2.704481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3041829) q[0];
sx q[0];
rz(-1.8619616) q[0];
sx q[0];
rz(0.090374723) q[0];
rz(-pi) q[1];
rz(-2.897981) q[2];
sx q[2];
rz(-1.1738452) q[2];
sx q[2];
rz(1.5543907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8902953) q[1];
sx q[1];
rz(-1.9922868) q[1];
sx q[1];
rz(-1.7226268) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0813619) q[3];
sx q[3];
rz(-2.0647991) q[3];
sx q[3];
rz(-1.4480818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1221216) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(1.9034889) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(-2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6195174) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(-2.1602901) q[0];
rz(0.13521067) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-0.19128004) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45368886) q[0];
sx q[0];
rz(-0.17246248) q[0];
sx q[0];
rz(1.0426636) q[0];
x q[1];
rz(0.27387597) q[2];
sx q[2];
rz(-0.855815) q[2];
sx q[2];
rz(0.82211923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8297255) q[1];
sx q[1];
rz(-2.3448181) q[1];
sx q[1];
rz(1.4273248) q[1];
x q[2];
rz(2.5883834) q[3];
sx q[3];
rz(-2.534453) q[3];
sx q[3];
rz(2.7943484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4613351) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(-0.7615532) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8028832) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(0.11511766) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(-2.1690878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2199729) q[0];
sx q[0];
rz(-0.12244206) q[0];
sx q[0];
rz(2.5570611) q[0];
x q[1];
rz(0.33072492) q[2];
sx q[2];
rz(-2.2933368) q[2];
sx q[2];
rz(-1.5028138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79341187) q[1];
sx q[1];
rz(-2.9017397) q[1];
sx q[1];
rz(2.1972149) q[1];
rz(0.036690849) q[3];
sx q[3];
rz(-0.9551691) q[3];
sx q[3];
rz(-0.43945593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82289034) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(-1.548432) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(2.1877066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71972972) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(1.4259889) q[0];
rz(-2.0772207) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(2.7672966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22898856) q[0];
sx q[0];
rz(-1.8870263) q[0];
sx q[0];
rz(2.7094748) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8315973) q[2];
sx q[2];
rz(-2.0506095) q[2];
sx q[2];
rz(-0.11670437) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4243851) q[1];
sx q[1];
rz(-1.8339515) q[1];
sx q[1];
rz(2.5724263) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2603287) q[3];
sx q[3];
rz(-2.1173819) q[3];
sx q[3];
rz(1.5342086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8950243) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(-0.22917497) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7763057) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(-2.2348485) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(1.3100756) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6543286) q[0];
sx q[0];
rz(-1.9582821) q[0];
sx q[0];
rz(-1.6261473) q[0];
rz(0.37808772) q[2];
sx q[2];
rz(-2.3921161) q[2];
sx q[2];
rz(2.7917002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8430427) q[1];
sx q[1];
rz(-1.7610465) q[1];
sx q[1];
rz(-2.802554) q[1];
rz(2.5861916) q[3];
sx q[3];
rz(-1.5941094) q[3];
sx q[3];
rz(-0.59613746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(1.6581992) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(-3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16335547) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(-0.50312463) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(-2.2917152) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06892878) q[0];
sx q[0];
rz(-1.7964296) q[0];
sx q[0];
rz(2.9805141) q[0];
x q[1];
rz(-0.88768994) q[2];
sx q[2];
rz(-1.6437093) q[2];
sx q[2];
rz(-2.2788252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0440158) q[1];
sx q[1];
rz(-0.8102639) q[1];
sx q[1];
rz(0.047659831) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0181846) q[3];
sx q[3];
rz(-1.3128237) q[3];
sx q[3];
rz(0.79515275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59051096) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(-1.3809416) q[2];
rz(2.3855709) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(-2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(-2.7365622) q[0];
rz(-0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.8639494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15690878) q[0];
sx q[0];
rz(-0.55159969) q[0];
sx q[0];
rz(-0.44996913) q[0];
rz(-pi) q[1];
rz(-0.69127609) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(0.28011766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9459878) q[1];
sx q[1];
rz(-2.0161649) q[1];
sx q[1];
rz(-3.1006378) q[1];
x q[2];
rz(0.38655917) q[3];
sx q[3];
rz(-0.87214008) q[3];
sx q[3];
rz(0.78702918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8032802) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(2.7837616) q[2];
rz(-1.4194277) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91838592) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(-0.11225587) q[0];
rz(2.2414801) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(2.9311438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4608599) q[0];
sx q[0];
rz(-2.0740777) q[0];
sx q[0];
rz(-1.3637278) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1591572) q[2];
sx q[2];
rz(-1.4201846) q[2];
sx q[2];
rz(2.6580236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.75281843) q[1];
sx q[1];
rz(-0.37467271) q[1];
sx q[1];
rz(2.545536) q[1];
rz(-pi) q[2];
rz(-0.25037346) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(1.9104513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9615053) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(1.5853184) q[2];
rz(1.2735584) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(0.20726985) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(-2.3251484) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(1.8893835) q[2];
sx q[2];
rz(-0.49124419) q[2];
sx q[2];
rz(-1.0791525) q[2];
rz(0.3420842) q[3];
sx q[3];
rz(-0.87089201) q[3];
sx q[3];
rz(1.0637829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
