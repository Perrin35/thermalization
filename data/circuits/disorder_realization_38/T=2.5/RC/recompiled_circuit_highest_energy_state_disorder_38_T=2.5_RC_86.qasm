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
rz(-2.3081245) q[0];
sx q[0];
rz(3.9552116) q[0];
sx q[0];
rz(7.8501346) q[0];
rz(-1.1235224) q[1];
sx q[1];
rz(-0.82668537) q[1];
sx q[1];
rz(-2.6348662) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0599466) q[0];
sx q[0];
rz(-1.5525408) q[0];
sx q[0];
rz(-0.55299802) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0000524) q[2];
sx q[2];
rz(-2.5876293) q[2];
sx q[2];
rz(0.20997722) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0953137) q[1];
sx q[1];
rz(-1.5161524) q[1];
sx q[1];
rz(-2.8846963) q[1];
rz(-0.039095239) q[3];
sx q[3];
rz(-2.2347898) q[3];
sx q[3];
rz(1.853626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.43164051) q[2];
sx q[2];
rz(-1.1974572) q[2];
sx q[2];
rz(2.6006827) q[2];
rz(-1.3276118) q[3];
sx q[3];
rz(-1.6874467) q[3];
sx q[3];
rz(0.66837627) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3427908) q[0];
sx q[0];
rz(-1.7690161) q[0];
sx q[0];
rz(2.0846682) q[0];
rz(0.42653719) q[1];
sx q[1];
rz(-1.6273472) q[1];
sx q[1];
rz(-3.1372517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89632242) q[0];
sx q[0];
rz(-1.9548804) q[0];
sx q[0];
rz(2.2784968) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7319381) q[2];
sx q[2];
rz(-0.97926295) q[2];
sx q[2];
rz(2.9681674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40758924) q[1];
sx q[1];
rz(-3.0619016) q[1];
sx q[1];
rz(1.7015639) q[1];
x q[2];
rz(2.6066296) q[3];
sx q[3];
rz(-1.6239979) q[3];
sx q[3];
rz(3.1411848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9428375) q[2];
sx q[2];
rz(-1.7386856) q[2];
sx q[2];
rz(-0.37503654) q[2];
rz(0.047920553) q[3];
sx q[3];
rz(-1.5791357) q[3];
sx q[3];
rz(-2.9722049) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7584943) q[0];
sx q[0];
rz(-2.5653745) q[0];
sx q[0];
rz(1.6297485) q[0];
rz(-1.0668628) q[1];
sx q[1];
rz(-1.8424415) q[1];
sx q[1];
rz(-2.3323434) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.209037) q[0];
sx q[0];
rz(-1.5309628) q[0];
sx q[0];
rz(-1.5341542) q[0];
rz(1.5009484) q[2];
sx q[2];
rz(-2.2874292) q[2];
sx q[2];
rz(1.6829862) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7748562) q[1];
sx q[1];
rz(-1.7637758) q[1];
sx q[1];
rz(-2.3425413) q[1];
rz(-1.567908) q[3];
sx q[3];
rz(-1.2089855) q[3];
sx q[3];
rz(-2.1199292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92144907) q[2];
sx q[2];
rz(-2.4704832) q[2];
sx q[2];
rz(-1.5044093) q[2];
rz(2.1680221) q[3];
sx q[3];
rz(-2.3928271) q[3];
sx q[3];
rz(-2.0481295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022738986) q[0];
sx q[0];
rz(-1.5819419) q[0];
sx q[0];
rz(-2.7795025) q[0];
rz(-1.7936961) q[1];
sx q[1];
rz(-1.7575512) q[1];
sx q[1];
rz(-1.0628343) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84447538) q[0];
sx q[0];
rz(-1.3985635) q[0];
sx q[0];
rz(-0.14396859) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9241829) q[2];
sx q[2];
rz(-2.1608781) q[2];
sx q[2];
rz(1.3979531) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3058051) q[1];
sx q[1];
rz(-1.9129681) q[1];
sx q[1];
rz(2.8306116) q[1];
x q[2];
rz(-1.3534668) q[3];
sx q[3];
rz(-1.3222851) q[3];
sx q[3];
rz(1.0729162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7245543) q[2];
sx q[2];
rz(-2.3596256) q[2];
sx q[2];
rz(0.057272591) q[2];
rz(0.96632424) q[3];
sx q[3];
rz(-0.94760197) q[3];
sx q[3];
rz(-1.9362601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22369497) q[0];
sx q[0];
rz(-0.88858336) q[0];
sx q[0];
rz(-2.5433871) q[0];
rz(-1.573645) q[1];
sx q[1];
rz(-1.6697829) q[1];
sx q[1];
rz(0.15984687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94598818) q[0];
sx q[0];
rz(-0.78778446) q[0];
sx q[0];
rz(0.24290084) q[0];
rz(-pi) q[1];
rz(2.8627928) q[2];
sx q[2];
rz(-2.9220916) q[2];
sx q[2];
rz(0.9448828) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5626403) q[1];
sx q[1];
rz(-0.75030164) q[1];
sx q[1];
rz(-2.5681488) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5261887) q[3];
sx q[3];
rz(-1.4945703) q[3];
sx q[3];
rz(-1.3520698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15888958) q[2];
sx q[2];
rz(-2.7003459) q[2];
sx q[2];
rz(-2.5510447) q[2];
rz(-0.96747893) q[3];
sx q[3];
rz(-1.5891821) q[3];
sx q[3];
rz(1.1837186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972873) q[0];
sx q[0];
rz(-0.66216457) q[0];
sx q[0];
rz(-1.039132) q[0];
rz(-0.93152535) q[1];
sx q[1];
rz(-2.4577699) q[1];
sx q[1];
rz(-1.1087803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88752969) q[0];
sx q[0];
rz(-1.0317804) q[0];
sx q[0];
rz(2.1655393) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80486561) q[2];
sx q[2];
rz(-2.1535843) q[2];
sx q[2];
rz(-1.031176) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1176374) q[1];
sx q[1];
rz(-1.9267356) q[1];
sx q[1];
rz(-0.61564095) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25082254) q[3];
sx q[3];
rz(-1.0366813) q[3];
sx q[3];
rz(-3.0297584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80591566) q[2];
sx q[2];
rz(-0.47199619) q[2];
sx q[2];
rz(-2.8373888) q[2];
rz(-0.16600569) q[3];
sx q[3];
rz(-1.7793572) q[3];
sx q[3];
rz(1.6148875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45749998) q[0];
sx q[0];
rz(-2.9170211) q[0];
sx q[0];
rz(-1.4373047) q[0];
rz(2.8999088) q[1];
sx q[1];
rz(-1.4983404) q[1];
sx q[1];
rz(-0.93745747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3110549) q[0];
sx q[0];
rz(-1.8094581) q[0];
sx q[0];
rz(-2.412459) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8514013) q[2];
sx q[2];
rz(-1.5414503) q[2];
sx q[2];
rz(1.9179232) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2087349) q[1];
sx q[1];
rz(-2.2003649) q[1];
sx q[1];
rz(2.1634794) q[1];
rz(-pi) q[2];
rz(-2.203905) q[3];
sx q[3];
rz(-2.4613639) q[3];
sx q[3];
rz(-0.89416789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9044056) q[2];
sx q[2];
rz(-2.0439456) q[2];
sx q[2];
rz(-2.912525) q[2];
rz(1.1711052) q[3];
sx q[3];
rz(-1.7541211) q[3];
sx q[3];
rz(0.68172014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1218629) q[0];
sx q[0];
rz(-0.64058146) q[0];
sx q[0];
rz(-1.4135452) q[0];
rz(1.2177263) q[1];
sx q[1];
rz(-1.7037337) q[1];
sx q[1];
rz(-0.53122836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41252688) q[0];
sx q[0];
rz(-2.7587037) q[0];
sx q[0];
rz(2.114243) q[0];
x q[1];
rz(-0.81905535) q[2];
sx q[2];
rz(-1.8764638) q[2];
sx q[2];
rz(1.3993939) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.574305) q[1];
sx q[1];
rz(-1.131408) q[1];
sx q[1];
rz(0.69171016) q[1];
x q[2];
rz(0.76721787) q[3];
sx q[3];
rz(-2.7306087) q[3];
sx q[3];
rz(2.1959636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56500498) q[2];
sx q[2];
rz(-1.7030623) q[2];
sx q[2];
rz(0.79565221) q[2];
rz(-2.3203881) q[3];
sx q[3];
rz(-2.932565) q[3];
sx q[3];
rz(2.8388099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32383701) q[0];
sx q[0];
rz(-2.9471286) q[0];
sx q[0];
rz(-1.6545464) q[0];
rz(1.5215634) q[1];
sx q[1];
rz(-1.227102) q[1];
sx q[1];
rz(-1.0618658) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.056039) q[0];
sx q[0];
rz(-2.0889707) q[0];
sx q[0];
rz(-2.0033625) q[0];
rz(-pi) q[1];
rz(0.050786917) q[2];
sx q[2];
rz(-1.2768695) q[2];
sx q[2];
rz(1.1370629) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1170845) q[1];
sx q[1];
rz(-1.7231795) q[1];
sx q[1];
rz(1.1445359) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2003593) q[3];
sx q[3];
rz(-1.9780156) q[3];
sx q[3];
rz(2.9349611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1294641) q[2];
sx q[2];
rz(-2.0120554) q[2];
sx q[2];
rz(-0.10936347) q[2];
rz(2.8520285) q[3];
sx q[3];
rz(-1.3907631) q[3];
sx q[3];
rz(0.65350986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43082032) q[0];
sx q[0];
rz(-1.1752952) q[0];
sx q[0];
rz(1.9868504) q[0];
rz(0.42354241) q[1];
sx q[1];
rz(-1.4332899) q[1];
sx q[1];
rz(0.67798859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2104312) q[0];
sx q[0];
rz(-0.95936459) q[0];
sx q[0];
rz(2.1290886) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7977436) q[2];
sx q[2];
rz(-1.316202) q[2];
sx q[2];
rz(2.7225021) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0667693) q[1];
sx q[1];
rz(-2.4500203) q[1];
sx q[1];
rz(-0.87705405) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5497694) q[3];
sx q[3];
rz(-1.4524959) q[3];
sx q[3];
rz(0.59676701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28422156) q[2];
sx q[2];
rz(-2.325433) q[2];
sx q[2];
rz(-0.59874272) q[2];
rz(1.5395928) q[3];
sx q[3];
rz(-1.2076104) q[3];
sx q[3];
rz(2.1046751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67958647) q[0];
sx q[0];
rz(-0.34272598) q[0];
sx q[0];
rz(1.7525679) q[0];
rz(-3.1287843) q[1];
sx q[1];
rz(-2.7638331) q[1];
sx q[1];
rz(1.4889181) q[1];
rz(-0.1107874) q[2];
sx q[2];
rz(-2.2089133) q[2];
sx q[2];
rz(-3.1089913) q[2];
rz(2.1297366) q[3];
sx q[3];
rz(-1.5717117) q[3];
sx q[3];
rz(-1.6757884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
