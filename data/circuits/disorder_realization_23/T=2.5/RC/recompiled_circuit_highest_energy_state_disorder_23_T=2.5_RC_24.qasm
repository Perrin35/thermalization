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
rz(-1.3462525) q[0];
sx q[0];
rz(-1.4486382) q[0];
sx q[0];
rz(1.1516655) q[0];
rz(1.0898074) q[1];
sx q[1];
rz(1.6648219) q[1];
sx q[1];
rz(9.2738168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1957832) q[0];
sx q[0];
rz(-1.8651143) q[0];
sx q[0];
rz(3.0867317) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8791228) q[2];
sx q[2];
rz(-2.8175857) q[2];
sx q[2];
rz(-0.25036795) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50569423) q[1];
sx q[1];
rz(-1.5855967) q[1];
sx q[1];
rz(3.1379051) q[1];
x q[2];
rz(0.096681194) q[3];
sx q[3];
rz(-1.5301609) q[3];
sx q[3];
rz(1.3717029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0207409) q[2];
sx q[2];
rz(-3.1181702) q[2];
sx q[2];
rz(2.6502996) q[2];
rz(1.322572) q[3];
sx q[3];
rz(-1.6072075) q[3];
sx q[3];
rz(1.3278495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.78111929) q[0];
sx q[0];
rz(-0.55447117) q[0];
sx q[0];
rz(-0.69960272) q[0];
rz(-1.5910925) q[1];
sx q[1];
rz(-0.49483776) q[1];
sx q[1];
rz(0.16986212) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3512974) q[0];
sx q[0];
rz(-1.5322432) q[0];
sx q[0];
rz(0.077705381) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3048068) q[2];
sx q[2];
rz(-1.6992603) q[2];
sx q[2];
rz(0.94513946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3292698) q[1];
sx q[1];
rz(-1.0877177) q[1];
sx q[1];
rz(0.4050576) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4109123) q[3];
sx q[3];
rz(-2.6794187) q[3];
sx q[3];
rz(1.9222207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8794609) q[2];
sx q[2];
rz(-0.36236557) q[2];
sx q[2];
rz(-1.7246838) q[2];
rz(-2.0587685) q[3];
sx q[3];
rz(-2.2542451) q[3];
sx q[3];
rz(0.82720238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5812434) q[0];
sx q[0];
rz(-2.2214948) q[0];
sx q[0];
rz(-1.5823407) q[0];
rz(0.70445591) q[1];
sx q[1];
rz(-1.1471986) q[1];
sx q[1];
rz(0.53680435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052633135) q[0];
sx q[0];
rz(-0.68872243) q[0];
sx q[0];
rz(3.0687276) q[0];
rz(-pi) q[1];
rz(3.051338) q[2];
sx q[2];
rz(-1.5343416) q[2];
sx q[2];
rz(-2.0608471) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0743064) q[1];
sx q[1];
rz(-0.18504772) q[1];
sx q[1];
rz(-0.61850278) q[1];
rz(1.8615253) q[3];
sx q[3];
rz(-0.71746263) q[3];
sx q[3];
rz(-0.91625253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.392936) q[2];
sx q[2];
rz(-1.3373249) q[2];
sx q[2];
rz(2.362747) q[2];
rz(-0.084176453) q[3];
sx q[3];
rz(-0.86936969) q[3];
sx q[3];
rz(1.4028153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8970784) q[0];
sx q[0];
rz(-2.2205413) q[0];
sx q[0];
rz(0.47065863) q[0];
rz(1.6828407) q[1];
sx q[1];
rz(-3.1367446) q[1];
sx q[1];
rz(-0.86476129) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1098561) q[0];
sx q[0];
rz(-1.512728) q[0];
sx q[0];
rz(0.51490291) q[0];
x q[1];
rz(2.607279) q[2];
sx q[2];
rz(-2.8861037) q[2];
sx q[2];
rz(1.8035672) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3833852) q[1];
sx q[1];
rz(-1.636519) q[1];
sx q[1];
rz(3.0905105) q[1];
x q[2];
rz(2.6506977) q[3];
sx q[3];
rz(-0.59661907) q[3];
sx q[3];
rz(2.9379549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.027815759) q[2];
sx q[2];
rz(-2.2769589) q[2];
sx q[2];
rz(1.9846385) q[2];
rz(-2.786934) q[3];
sx q[3];
rz(-1.2361453) q[3];
sx q[3];
rz(-0.11053301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81195867) q[0];
sx q[0];
rz(-0.93557731) q[0];
sx q[0];
rz(-0.53363609) q[0];
rz(1.927902) q[1];
sx q[1];
rz(-0.023093725) q[1];
sx q[1];
rz(1.1812706) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1071067) q[0];
sx q[0];
rz(-1.6167322) q[0];
sx q[0];
rz(-0.1631384) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1540765) q[2];
sx q[2];
rz(-0.96229711) q[2];
sx q[2];
rz(1.6241983) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5983728) q[1];
sx q[1];
rz(-0.66987619) q[1];
sx q[1];
rz(1.8240488) q[1];
rz(-pi) q[2];
rz(1.817068) q[3];
sx q[3];
rz(-2.6239873) q[3];
sx q[3];
rz(0.31806614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10826762) q[2];
sx q[2];
rz(-0.51563087) q[2];
sx q[2];
rz(0.92833129) q[2];
rz(-0.27836529) q[3];
sx q[3];
rz(-1.8230349) q[3];
sx q[3];
rz(3.0729455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5334897) q[0];
sx q[0];
rz(-2.6187496) q[0];
sx q[0];
rz(-2.9485517) q[0];
rz(-0.67735425) q[1];
sx q[1];
rz(-0.028379863) q[1];
sx q[1];
rz(-1.5100286) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91514912) q[0];
sx q[0];
rz(-2.9993176) q[0];
sx q[0];
rz(2.9998542) q[0];
rz(-pi) q[1];
rz(-1.6608428) q[2];
sx q[2];
rz(-0.46381942) q[2];
sx q[2];
rz(-0.12741379) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7379166) q[1];
sx q[1];
rz(-2.3737646) q[1];
sx q[1];
rz(-1.9270792) q[1];
rz(-pi) q[2];
rz(-2.5872748) q[3];
sx q[3];
rz(-0.80762562) q[3];
sx q[3];
rz(0.19632158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3626927) q[2];
sx q[2];
rz(-0.71234667) q[2];
sx q[2];
rz(-2.1544797) q[2];
rz(1.0846064) q[3];
sx q[3];
rz(-2.5955213) q[3];
sx q[3];
rz(2.5206821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5622332) q[0];
sx q[0];
rz(-2.1491304) q[0];
sx q[0];
rz(1.5935422) q[0];
rz(-0.69001895) q[1];
sx q[1];
rz(-0.023612173) q[1];
sx q[1];
rz(1.7049559) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6520335) q[0];
sx q[0];
rz(-1.0512182) q[0];
sx q[0];
rz(-1.621031) q[0];
rz(-pi) q[1];
x q[1];
rz(0.026937303) q[2];
sx q[2];
rz(-0.68490324) q[2];
sx q[2];
rz(1.2482582) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44013906) q[1];
sx q[1];
rz(-1.6303282) q[1];
sx q[1];
rz(-2.1131705) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5593188) q[3];
sx q[3];
rz(-0.52460656) q[3];
sx q[3];
rz(-0.23422262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37316608) q[2];
sx q[2];
rz(-1.2370141) q[2];
sx q[2];
rz(2.5567059) q[2];
rz(-1.5960426) q[3];
sx q[3];
rz(-1.7187748) q[3];
sx q[3];
rz(-2.2271631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5717413) q[0];
sx q[0];
rz(-2.2241156) q[0];
sx q[0];
rz(1.5745987) q[0];
rz(0.51708108) q[1];
sx q[1];
rz(-2.2061429) q[1];
sx q[1];
rz(1.8749974) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168473) q[0];
sx q[0];
rz(-3.1355097) q[0];
sx q[0];
rz(-2.541126) q[0];
x q[1];
rz(0.80516239) q[2];
sx q[2];
rz(-1.9765719) q[2];
sx q[2];
rz(-2.9540887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4609485) q[1];
sx q[1];
rz(-1.0310182) q[1];
sx q[1];
rz(-2.7311027) q[1];
rz(-pi) q[2];
rz(1.6932159) q[3];
sx q[3];
rz(-1.4493692) q[3];
sx q[3];
rz(-1.9521015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7685585) q[2];
sx q[2];
rz(-0.80058241) q[2];
sx q[2];
rz(1.2726834) q[2];
rz(-0.4396762) q[3];
sx q[3];
rz(-1.9821854) q[3];
sx q[3];
rz(3.0769949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8398297) q[0];
sx q[0];
rz(-3.1253212) q[0];
sx q[0];
rz(-2.8453258) q[0];
rz(-0.06427327) q[1];
sx q[1];
rz(-1.6769033) q[1];
sx q[1];
rz(1.6305264) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13416269) q[0];
sx q[0];
rz(-1.2399956) q[0];
sx q[0];
rz(-0.13749977) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4687541) q[2];
sx q[2];
rz(-1.7305859) q[2];
sx q[2];
rz(0.21110134) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8478197) q[1];
sx q[1];
rz(-2.7917531) q[1];
sx q[1];
rz(1.185525) q[1];
rz(-1.4496719) q[3];
sx q[3];
rz(-1.559404) q[3];
sx q[3];
rz(-0.73339547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5955547) q[2];
sx q[2];
rz(-1.7067355) q[2];
sx q[2];
rz(2.0175374) q[2];
rz(1.837364) q[3];
sx q[3];
rz(-0.29574695) q[3];
sx q[3];
rz(-0.058163253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.27125636) q[0];
sx q[0];
rz(-2.3503292) q[0];
sx q[0];
rz(-1.9747718) q[0];
rz(-1.600949) q[1];
sx q[1];
rz(-0.29778844) q[1];
sx q[1];
rz(-1.8150394) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2873424) q[0];
sx q[0];
rz(-2.4357987) q[0];
sx q[0];
rz(2.7588034) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9341078) q[2];
sx q[2];
rz(-1.9544387) q[2];
sx q[2];
rz(0.07019474) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1455166) q[1];
sx q[1];
rz(-2.4156486) q[1];
sx q[1];
rz(-0.32246253) q[1];
x q[2];
rz(0.097685944) q[3];
sx q[3];
rz(-1.9703416) q[3];
sx q[3];
rz(2.0942502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88663179) q[2];
sx q[2];
rz(-0.16356629) q[2];
sx q[2];
rz(-1.4476042) q[2];
rz(3.0149031) q[3];
sx q[3];
rz(-1.6619752) q[3];
sx q[3];
rz(-1.1160342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3348677) q[0];
sx q[0];
rz(-1.3074449) q[0];
sx q[0];
rz(1.773651) q[0];
rz(-1.5927636) q[1];
sx q[1];
rz(-0.82850414) q[1];
sx q[1];
rz(-2.9864476) q[1];
rz(-0.44687475) q[2];
sx q[2];
rz(-1.7394273) q[2];
sx q[2];
rz(1.5706825) q[2];
rz(2.9709076) q[3];
sx q[3];
rz(-2.1483834) q[3];
sx q[3];
rz(-2.7494062) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
