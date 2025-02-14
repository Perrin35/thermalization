OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(4.4409039) q[0];
sx q[0];
rz(9.470603) q[0];
rz(-1.8419645) q[1];
sx q[1];
rz(-1.3388495) q[1];
sx q[1];
rz(-1.2531228) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1880671) q[0];
sx q[0];
rz(-1.8154411) q[0];
sx q[0];
rz(1.1172386) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21109493) q[2];
sx q[2];
rz(-2.935754) q[2];
sx q[2];
rz(-1.7072276) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7980849) q[1];
sx q[1];
rz(-1.6899334) q[1];
sx q[1];
rz(2.3878674) q[1];
rz(-pi) q[2];
rz(-1.5398938) q[3];
sx q[3];
rz(-0.12750827) q[3];
sx q[3];
rz(2.7650583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.32156285) q[2];
sx q[2];
rz(-1.7367312) q[2];
sx q[2];
rz(0.31537867) q[2];
rz(-0.086221181) q[3];
sx q[3];
rz(-3.0486221) q[3];
sx q[3];
rz(2.2975547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.740199) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(-0.33541086) q[0];
rz(2.7566578) q[1];
sx q[1];
rz(-0.79610577) q[1];
sx q[1];
rz(1.2892494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6008284) q[0];
sx q[0];
rz(-2.7473831) q[0];
sx q[0];
rz(2.9790131) q[0];
x q[1];
rz(-2.2122967) q[2];
sx q[2];
rz(-1.6469064) q[2];
sx q[2];
rz(0.92744213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1114167) q[1];
sx q[1];
rz(-1.1892288) q[1];
sx q[1];
rz(-1.2658582) q[1];
rz(-0.13522526) q[3];
sx q[3];
rz(-1.8117419) q[3];
sx q[3];
rz(1.2775482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5424767) q[2];
sx q[2];
rz(-1.7581538) q[2];
sx q[2];
rz(1.4625589) q[2];
rz(2.2595432) q[3];
sx q[3];
rz(-2.4817011) q[3];
sx q[3];
rz(1.6409469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3682692) q[0];
sx q[0];
rz(-2.6550846) q[0];
sx q[0];
rz(2.4988556) q[0];
rz(2.2679988) q[1];
sx q[1];
rz(-1.3106376) q[1];
sx q[1];
rz(-2.1547623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8514252) q[0];
sx q[0];
rz(-2.254524) q[0];
sx q[0];
rz(-0.19398035) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6880694) q[2];
sx q[2];
rz(-2.6542695) q[2];
sx q[2];
rz(0.27299689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3915186) q[1];
sx q[1];
rz(-1.0315391) q[1];
sx q[1];
rz(-2.2184994) q[1];
x q[2];
rz(0.049792265) q[3];
sx q[3];
rz(-1.9450499) q[3];
sx q[3];
rz(-0.58993898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4387536) q[2];
sx q[2];
rz(-3.1162016) q[2];
sx q[2];
rz(-2.0849483) q[2];
rz(1.3422525) q[3];
sx q[3];
rz(-1.6868351) q[3];
sx q[3];
rz(2.2630283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037647) q[0];
sx q[0];
rz(-2.8570638) q[0];
sx q[0];
rz(-1.1430662) q[0];
rz(-0.68483886) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(0.24982223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78501084) q[0];
sx q[0];
rz(-2.2306721) q[0];
sx q[0];
rz(2.5500245) q[0];
rz(-pi) q[1];
rz(1.5748137) q[2];
sx q[2];
rz(-2.145014) q[2];
sx q[2];
rz(-1.9030273) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64678363) q[1];
sx q[1];
rz(-0.80776359) q[1];
sx q[1];
rz(2.450472) q[1];
x q[2];
rz(1.8127727) q[3];
sx q[3];
rz(-1.5692838) q[3];
sx q[3];
rz(0.63336271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7217241) q[2];
sx q[2];
rz(-1.7896174) q[2];
sx q[2];
rz(-1.8458337) q[2];
rz(1.7660247) q[3];
sx q[3];
rz(-1.4294759) q[3];
sx q[3];
rz(-3.0425369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38882035) q[0];
sx q[0];
rz(-1.7701912) q[0];
sx q[0];
rz(0.8465299) q[0];
rz(-1.1835774) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(-1.2394946) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13279937) q[0];
sx q[0];
rz(-1.4504878) q[0];
sx q[0];
rz(1.3161206) q[0];
rz(2.7991512) q[2];
sx q[2];
rz(-1.3991809) q[2];
sx q[2];
rz(-1.1521074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2927276) q[1];
sx q[1];
rz(-1.5827109) q[1];
sx q[1];
rz(1.5148276) q[1];
rz(3.0817075) q[3];
sx q[3];
rz(-2.1671173) q[3];
sx q[3];
rz(-2.8879762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6401297) q[2];
sx q[2];
rz(-1.5093466) q[2];
sx q[2];
rz(-0.35422361) q[2];
rz(-0.7192449) q[3];
sx q[3];
rz(-2.5651599) q[3];
sx q[3];
rz(2.332212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7648014) q[0];
sx q[0];
rz(-0.23623315) q[0];
sx q[0];
rz(-0.34673044) q[0];
rz(-0.72225839) q[1];
sx q[1];
rz(-1.5303231) q[1];
sx q[1];
rz(2.9647656) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2472025) q[0];
sx q[0];
rz(-0.48485928) q[0];
sx q[0];
rz(1.2653744) q[0];
rz(2.0697303) q[2];
sx q[2];
rz(-1.6526573) q[2];
sx q[2];
rz(-2.6112542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7737501) q[1];
sx q[1];
rz(-2.9108139) q[1];
sx q[1];
rz(2.6683067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7309639) q[3];
sx q[3];
rz(-2.5729542) q[3];
sx q[3];
rz(-1.4071495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5727545) q[2];
sx q[2];
rz(-0.73358959) q[2];
sx q[2];
rz(2.0983762) q[2];
rz(-2.9853232) q[3];
sx q[3];
rz(-1.0633435) q[3];
sx q[3];
rz(-0.094303057) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25065502) q[0];
sx q[0];
rz(-1.0303048) q[0];
sx q[0];
rz(1.913273) q[0];
rz(0.43371513) q[1];
sx q[1];
rz(-0.80077306) q[1];
sx q[1];
rz(2.0965651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21711536) q[0];
sx q[0];
rz(-1.7733367) q[0];
sx q[0];
rz(3.0826352) q[0];
rz(2.074998) q[2];
sx q[2];
rz(-0.81461755) q[2];
sx q[2];
rz(0.93570081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8974287) q[1];
sx q[1];
rz(-2.5153195) q[1];
sx q[1];
rz(1.3757964) q[1];
rz(-0.52900903) q[3];
sx q[3];
rz(-2.2914791) q[3];
sx q[3];
rz(2.1095654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9016483) q[2];
sx q[2];
rz(-1.0239235) q[2];
sx q[2];
rz(0.18590064) q[2];
rz(1.2402041) q[3];
sx q[3];
rz(-1.5408206) q[3];
sx q[3];
rz(-2.127229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1321201) q[0];
sx q[0];
rz(-1.2484231) q[0];
sx q[0];
rz(-0.15952071) q[0];
rz(-2.3347143) q[1];
sx q[1];
rz(-1.2346377) q[1];
sx q[1];
rz(3.1330915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91220605) q[0];
sx q[0];
rz(-1.1649781) q[0];
sx q[0];
rz(-0.1677558) q[0];
rz(1.6775756) q[2];
sx q[2];
rz(-2.3651383) q[2];
sx q[2];
rz(-1.8242893) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8574243) q[1];
sx q[1];
rz(-1.4393974) q[1];
sx q[1];
rz(-0.22716503) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0866585) q[3];
sx q[3];
rz(-1.7268506) q[3];
sx q[3];
rz(-2.597996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8622417) q[2];
sx q[2];
rz(-1.3631577) q[2];
sx q[2];
rz(-0.19374338) q[2];
rz(1.6019609) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(2.0120473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2368471) q[0];
sx q[0];
rz(-1.756825) q[0];
sx q[0];
rz(0.73614502) q[0];
rz(-0.88788095) q[1];
sx q[1];
rz(-0.45279756) q[1];
sx q[1];
rz(-3.091541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5745478) q[0];
sx q[0];
rz(-2.7549681) q[0];
sx q[0];
rz(2.6599314) q[0];
rz(-2.8456837) q[2];
sx q[2];
rz(-0.32797932) q[2];
sx q[2];
rz(0.28959549) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.012573012) q[1];
sx q[1];
rz(-0.37380773) q[1];
sx q[1];
rz(0.5859979) q[1];
rz(-pi) q[2];
rz(0.56882755) q[3];
sx q[3];
rz(-2.1095368) q[3];
sx q[3];
rz(0.26876581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8686409) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(0.0090553332) q[2];
rz(-2.791259) q[3];
sx q[3];
rz(-2.83559) q[3];
sx q[3];
rz(-0.30456021) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41651273) q[0];
sx q[0];
rz(-2.0820936) q[0];
sx q[0];
rz(-0.9730202) q[0];
rz(1.579772) q[1];
sx q[1];
rz(-2.2353954) q[1];
sx q[1];
rz(-2.3320893) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1081297) q[0];
sx q[0];
rz(-1.9282883) q[0];
sx q[0];
rz(-0.78127677) q[0];
x q[1];
rz(-2.948582) q[2];
sx q[2];
rz(-1.7959611) q[2];
sx q[2];
rz(-1.7042773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5132147) q[1];
sx q[1];
rz(-0.49114463) q[1];
sx q[1];
rz(-1.9188736) q[1];
rz(-pi) q[2];
rz(-2.2274251) q[3];
sx q[3];
rz(-0.69147136) q[3];
sx q[3];
rz(-1.4007614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1805264) q[2];
sx q[2];
rz(-2.3193391) q[2];
sx q[2];
rz(2.6194438) q[2];
rz(2.7850049) q[3];
sx q[3];
rz(-0.57727376) q[3];
sx q[3];
rz(3.0780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4840354) q[0];
sx q[0];
rz(-0.35150305) q[0];
sx q[0];
rz(1.2765314) q[0];
rz(-0.36318489) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(0.77728669) q[2];
sx q[2];
rz(-1.9718134) q[2];
sx q[2];
rz(2.8169463) q[2];
rz(-1.894886) q[3];
sx q[3];
rz(-2.1264599) q[3];
sx q[3];
rz(2.651941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
