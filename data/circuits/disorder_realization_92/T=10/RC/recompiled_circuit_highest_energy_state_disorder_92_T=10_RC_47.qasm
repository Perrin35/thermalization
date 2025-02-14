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
rz(-3.0109316) q[0];
sx q[0];
rz(-1.8008512) q[0];
sx q[0];
rz(-1.2641579) q[0];
rz(-2.2364927) q[1];
sx q[1];
rz(-2.0533419) q[1];
sx q[1];
rz(-0.01297125) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.599975) q[0];
sx q[0];
rz(-2.0429039) q[0];
sx q[0];
rz(-1.3946482) q[0];
rz(-pi) q[1];
rz(1.4616724) q[2];
sx q[2];
rz(-1.4089917) q[2];
sx q[2];
rz(-1.2798123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2544059) q[1];
sx q[1];
rz(-1.2206843) q[1];
sx q[1];
rz(1.2302047) q[1];
x q[2];
rz(-1.3566586) q[3];
sx q[3];
rz(-0.99839568) q[3];
sx q[3];
rz(2.1635087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26244792) q[2];
sx q[2];
rz(-1.6877354) q[2];
sx q[2];
rz(-0.095414735) q[2];
rz(-0.48424193) q[3];
sx q[3];
rz(-0.80317322) q[3];
sx q[3];
rz(1.6835083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478304) q[0];
sx q[0];
rz(-1.7050803) q[0];
sx q[0];
rz(2.4875212) q[0];
rz(1.7895128) q[1];
sx q[1];
rz(-0.65579954) q[1];
sx q[1];
rz(0.10993122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9655361) q[0];
sx q[0];
rz(-1.5829464) q[0];
sx q[0];
rz(1.5506844) q[0];
rz(-0.45808001) q[2];
sx q[2];
rz(-0.9914248) q[2];
sx q[2];
rz(1.9922076) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.667404) q[1];
sx q[1];
rz(-1.1193573) q[1];
sx q[1];
rz(-1.3787446) q[1];
x q[2];
rz(2.7450493) q[3];
sx q[3];
rz(-1.8064883) q[3];
sx q[3];
rz(-1.1514219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.316651) q[2];
sx q[2];
rz(-1.2085088) q[2];
sx q[2];
rz(1.6744772) q[2];
rz(2.1064827) q[3];
sx q[3];
rz(-2.3605774) q[3];
sx q[3];
rz(1.8234113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255945) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(-0.79063928) q[0];
rz(0.74360338) q[1];
sx q[1];
rz(-1.598282) q[1];
sx q[1];
rz(-0.57949439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6181439) q[0];
sx q[0];
rz(-2.14144) q[0];
sx q[0];
rz(-2.4916925) q[0];
x q[1];
rz(-1.2762464) q[2];
sx q[2];
rz(-1.8229457) q[2];
sx q[2];
rz(-1.1974826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6404785) q[1];
sx q[1];
rz(-1.6081297) q[1];
sx q[1];
rz(-1.3126612) q[1];
rz(-pi) q[2];
rz(2.0602588) q[3];
sx q[3];
rz(-2.0474985) q[3];
sx q[3];
rz(2.7923194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6774595) q[2];
sx q[2];
rz(-0.83659283) q[2];
sx q[2];
rz(0.38132384) q[2];
rz(2.2903806) q[3];
sx q[3];
rz(-2.2150453) q[3];
sx q[3];
rz(-0.73494953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6094991) q[0];
sx q[0];
rz(-2.5272326) q[0];
sx q[0];
rz(2.7771948) q[0];
rz(0.20026194) q[1];
sx q[1];
rz(-1.4118782) q[1];
sx q[1];
rz(3.0416378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0138977) q[0];
sx q[0];
rz(-1.5504) q[0];
sx q[0];
rz(-1.4958032) q[0];
x q[1];
rz(0.1047524) q[2];
sx q[2];
rz(-1.9974553) q[2];
sx q[2];
rz(1.1327281) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2435088) q[1];
sx q[1];
rz(-1.5931411) q[1];
sx q[1];
rz(-0.0925272) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95619802) q[3];
sx q[3];
rz(-2.0609988) q[3];
sx q[3];
rz(-2.4407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0692856) q[2];
sx q[2];
rz(-1.0682718) q[2];
sx q[2];
rz(-3.0266673) q[2];
rz(-1.140444) q[3];
sx q[3];
rz(-1.2830257) q[3];
sx q[3];
rz(2.1663402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045778) q[0];
sx q[0];
rz(-0.95252043) q[0];
sx q[0];
rz(0.17247795) q[0];
rz(2.1306439) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(-0.15484658) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2456937) q[0];
sx q[0];
rz(-1.2749199) q[0];
sx q[0];
rz(2.9135743) q[0];
x q[1];
rz(-1.2482951) q[2];
sx q[2];
rz(-0.66993827) q[2];
sx q[2];
rz(0.33316406) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6754969) q[1];
sx q[1];
rz(-1.7148397) q[1];
sx q[1];
rz(-2.7269468) q[1];
x q[2];
rz(1.6825292) q[3];
sx q[3];
rz(-1.3936685) q[3];
sx q[3];
rz(2.2551409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0871206) q[2];
sx q[2];
rz(-2.8643769) q[2];
sx q[2];
rz(0.37230125) q[2];
rz(2.7436658) q[3];
sx q[3];
rz(-1.9979265) q[3];
sx q[3];
rz(-3.1411689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18846866) q[0];
sx q[0];
rz(-2.5690881) q[0];
sx q[0];
rz(1.5931801) q[0];
rz(0.36987034) q[1];
sx q[1];
rz(-1.3984171) q[1];
sx q[1];
rz(-0.63873783) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34015481) q[0];
sx q[0];
rz(-1.7246913) q[0];
sx q[0];
rz(1.3305386) q[0];
x q[1];
rz(0.20301007) q[2];
sx q[2];
rz(-1.0542717) q[2];
sx q[2];
rz(0.85337108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46203278) q[1];
sx q[1];
rz(-0.61374329) q[1];
sx q[1];
rz(-0.8590974) q[1];
x q[2];
rz(-0.79213284) q[3];
sx q[3];
rz(-1.9895305) q[3];
sx q[3];
rz(1.2611539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0729735) q[2];
sx q[2];
rz(-1.0516473) q[2];
sx q[2];
rz(-0.2529141) q[2];
rz(-0.84826338) q[3];
sx q[3];
rz(-1.6631283) q[3];
sx q[3];
rz(0.084913582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0406168) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(-2.9926391) q[0];
rz(-1.8303998) q[1];
sx q[1];
rz(-1.4102178) q[1];
sx q[1];
rz(-0.054100903) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9210756) q[0];
sx q[0];
rz(-2.069325) q[0];
sx q[0];
rz(-1.6706628) q[0];
rz(-pi) q[1];
rz(-3.0981423) q[2];
sx q[2];
rz(-0.80294007) q[2];
sx q[2];
rz(-2.174365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5230519) q[1];
sx q[1];
rz(-2.4344011) q[1];
sx q[1];
rz(-2.9721353) q[1];
x q[2];
rz(3.1182958) q[3];
sx q[3];
rz(-1.0956278) q[3];
sx q[3];
rz(3.0123229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3397303) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(-2.1745963) q[2];
rz(0.70080924) q[3];
sx q[3];
rz(-2.1613224) q[3];
sx q[3];
rz(2.7907659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
rz(0.12857777) q[0];
sx q[0];
rz(-1.9676493) q[0];
sx q[0];
rz(1.5482192) q[0];
rz(0.37823996) q[1];
sx q[1];
rz(-2.389237) q[1];
sx q[1];
rz(-0.38984782) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8075832) q[0];
sx q[0];
rz(-2.5981059) q[0];
sx q[0];
rz(-0.38477151) q[0];
rz(-pi) q[1];
rz(1.5178096) q[2];
sx q[2];
rz(-0.30445004) q[2];
sx q[2];
rz(-1.9805769) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88943849) q[1];
sx q[1];
rz(-0.61823003) q[1];
sx q[1];
rz(0.63207027) q[1];
rz(-pi) q[2];
rz(-2.8865783) q[3];
sx q[3];
rz(-1.7863635) q[3];
sx q[3];
rz(-0.6620342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0491911) q[2];
sx q[2];
rz(-2.0202049) q[2];
sx q[2];
rz(-1.8828877) q[2];
rz(0.24426584) q[3];
sx q[3];
rz(-1.786307) q[3];
sx q[3];
rz(-2.145483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8213537) q[0];
sx q[0];
rz(-1.2293674) q[0];
sx q[0];
rz(1.7852596) q[0];
rz(-2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(-2.70539) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2552196) q[0];
sx q[0];
rz(-2.5804735) q[0];
sx q[0];
rz(2.5291689) q[0];
x q[1];
rz(-0.8955375) q[2];
sx q[2];
rz(-2.6288599) q[2];
sx q[2];
rz(0.89287478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6375745) q[1];
sx q[1];
rz(-2.3406041) q[1];
sx q[1];
rz(1.8941419) q[1];
rz(-0.14472503) q[3];
sx q[3];
rz(-0.44071482) q[3];
sx q[3];
rz(-3.0851229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1205552) q[2];
sx q[2];
rz(-1.977481) q[2];
sx q[2];
rz(-2.5980921) q[2];
rz(-0.049526878) q[3];
sx q[3];
rz(-1.6560358) q[3];
sx q[3];
rz(1.4878976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5307584) q[0];
sx q[0];
rz(-3.0093091) q[0];
sx q[0];
rz(-3.0623867) q[0];
rz(2.7414956) q[1];
sx q[1];
rz(-2.2875417) q[1];
sx q[1];
rz(1.0265464) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8456789) q[0];
sx q[0];
rz(-2.2845553) q[0];
sx q[0];
rz(0.31166021) q[0];
rz(-1.4756938) q[2];
sx q[2];
rz(-1.9355023) q[2];
sx q[2];
rz(-2.022418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0821918) q[1];
sx q[1];
rz(-2.804356) q[1];
sx q[1];
rz(0.2880917) q[1];
x q[2];
rz(-0.35086326) q[3];
sx q[3];
rz(-1.6454344) q[3];
sx q[3];
rz(0.78001744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2058699) q[2];
sx q[2];
rz(-1.3047855) q[2];
sx q[2];
rz(1.4523466) q[2];
rz(0.82990372) q[3];
sx q[3];
rz(-2.2645576) q[3];
sx q[3];
rz(-0.95480603) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6185388) q[0];
sx q[0];
rz(-1.1445615) q[0];
sx q[0];
rz(-1.6339697) q[0];
rz(-1.2670831) q[1];
sx q[1];
rz(-1.0825842) q[1];
sx q[1];
rz(-2.4172197) q[1];
rz(0.53955033) q[2];
sx q[2];
rz(-1.0536516) q[2];
sx q[2];
rz(-0.45894844) q[2];
rz(-1.2290365) q[3];
sx q[3];
rz(-2.1038341) q[3];
sx q[3];
rz(-2.9084222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
