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
rz(1.4003657) q[0];
sx q[0];
rz(-0.23973149) q[0];
sx q[0];
rz(0.32455197) q[0];
rz(5.3311081) q[1];
sx q[1];
rz(6.0573112) q[1];
sx q[1];
rz(4.4499302) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1150779) q[0];
sx q[0];
rz(-2.240772) q[0];
sx q[0];
rz(1.3967394) q[0];
rz(-1.4412854) q[2];
sx q[2];
rz(-2.0515832) q[2];
sx q[2];
rz(-0.12809424) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0066915) q[1];
sx q[1];
rz(-2.6009702) q[1];
sx q[1];
rz(2.0903793) q[1];
rz(-pi) q[2];
rz(-2.9746858) q[3];
sx q[3];
rz(-2.4354773) q[3];
sx q[3];
rz(-2.0049068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41210458) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(-0.35120249) q[2];
rz(1.2900194) q[3];
sx q[3];
rz(-2.0100287) q[3];
sx q[3];
rz(1.2235519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38607645) q[0];
sx q[0];
rz(-1.9367243) q[0];
sx q[0];
rz(2.2112041) q[0];
rz(1.3049841) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(-2.5321541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4131282) q[0];
sx q[0];
rz(-0.43914686) q[0];
sx q[0];
rz(-1.8154316) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2728884) q[2];
sx q[2];
rz(-1.6570083) q[2];
sx q[2];
rz(2.6484368) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6516782) q[1];
sx q[1];
rz(-0.53333542) q[1];
sx q[1];
rz(-0.79483219) q[1];
rz(-pi) q[2];
rz(-2.9744963) q[3];
sx q[3];
rz(-1.0067938) q[3];
sx q[3];
rz(-0.43844863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0448407) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(0.074782221) q[2];
rz(-2.6212202) q[3];
sx q[3];
rz(-2.4986391) q[3];
sx q[3];
rz(0.61409942) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5388913) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-0.30211788) q[0];
rz(2.865454) q[1];
sx q[1];
rz(-1.8243676) q[1];
sx q[1];
rz(1.709323) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30423588) q[0];
sx q[0];
rz(-0.68528995) q[0];
sx q[0];
rz(-2.1987112) q[0];
rz(-pi) q[1];
rz(1.309607) q[2];
sx q[2];
rz(-0.99203324) q[2];
sx q[2];
rz(2.8566993) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.182309) q[1];
sx q[1];
rz(-1.4136836) q[1];
sx q[1];
rz(-1.0090488) q[1];
rz(-3.0705419) q[3];
sx q[3];
rz(-2.0356352) q[3];
sx q[3];
rz(-2.0185061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8351195) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(1.845537) q[2];
rz(-2.6895788) q[3];
sx q[3];
rz(-2.0343503) q[3];
sx q[3];
rz(-0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6119659) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(1.8091328) q[0];
rz(-2.0924856) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(-1.5886935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6352672) q[0];
sx q[0];
rz(-0.46231368) q[0];
sx q[0];
rz(0.72700951) q[0];
x q[1];
rz(2.0715782) q[2];
sx q[2];
rz(-1.5827279) q[2];
sx q[2];
rz(-2.3119761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59999787) q[1];
sx q[1];
rz(-1.7686971) q[1];
sx q[1];
rz(2.2210414) q[1];
rz(-2.1721607) q[3];
sx q[3];
rz(-2.3323389) q[3];
sx q[3];
rz(-1.7084029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66854746) q[2];
sx q[2];
rz(-2.1210402) q[2];
sx q[2];
rz(-0.2505396) q[2];
rz(-0.9969095) q[3];
sx q[3];
rz(-2.0018115) q[3];
sx q[3];
rz(-2.7610049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8119891) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(-1.3478152) q[0];
rz(0.40428058) q[1];
sx q[1];
rz(-0.75899044) q[1];
sx q[1];
rz(2.0124729) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383352) q[0];
sx q[0];
rz(-1.9651439) q[0];
sx q[0];
rz(0.82259891) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90520827) q[2];
sx q[2];
rz(-1.357175) q[2];
sx q[2];
rz(2.7366432) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3394062) q[1];
sx q[1];
rz(-1.4660343) q[1];
sx q[1];
rz(2.2441007) q[1];
rz(-pi) q[2];
rz(-2.8267838) q[3];
sx q[3];
rz(-1.6757351) q[3];
sx q[3];
rz(-2.2887857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3471442) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(1.8974737) q[2];
rz(-1.0602903) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4853972) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(0.42718497) q[0];
rz(-3.1071013) q[1];
sx q[1];
rz(-1.5723615) q[1];
sx q[1];
rz(0.010206612) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611864) q[0];
sx q[0];
rz(-1.5725333) q[0];
sx q[0];
rz(0.0023567452) q[0];
rz(-pi) q[1];
rz(0.16049196) q[2];
sx q[2];
rz(-0.94466034) q[2];
sx q[2];
rz(2.6662835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31406826) q[1];
sx q[1];
rz(-1.5064824) q[1];
sx q[1];
rz(-2.2576648) q[1];
x q[2];
rz(2.4236492) q[3];
sx q[3];
rz(-0.69980757) q[3];
sx q[3];
rz(-1.0216624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61937845) q[2];
sx q[2];
rz(-2.7284315) q[2];
sx q[2];
rz(-0.79356066) q[2];
rz(0.86269745) q[3];
sx q[3];
rz(-1.5669275) q[3];
sx q[3];
rz(2.4376552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4392387) q[0];
sx q[0];
rz(-2.2965501) q[0];
sx q[0];
rz(2.0080361) q[0];
rz(-1.0776862) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(0.30212197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0724807) q[0];
sx q[0];
rz(-1.9078507) q[0];
sx q[0];
rz(1.891579) q[0];
rz(2.0321192) q[2];
sx q[2];
rz(-1.0240842) q[2];
sx q[2];
rz(-0.49671587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6985642) q[1];
sx q[1];
rz(-2.5686567) q[1];
sx q[1];
rz(-0.40614508) q[1];
rz(-pi) q[2];
rz(-0.68639836) q[3];
sx q[3];
rz(-1.0960311) q[3];
sx q[3];
rz(1.7976185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2332396) q[2];
sx q[2];
rz(-2.2519604) q[2];
sx q[2];
rz(1.7944149) q[2];
rz(1.8917278) q[3];
sx q[3];
rz(-1.1976539) q[3];
sx q[3];
rz(-3.0344322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5477448) q[0];
sx q[0];
rz(-2.7084454) q[0];
sx q[0];
rz(3.0331842) q[0];
rz(1.0477061) q[1];
sx q[1];
rz(-2.3240418) q[1];
sx q[1];
rz(-1.0188867) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86025809) q[0];
sx q[0];
rz(-1.0509914) q[0];
sx q[0];
rz(-0.61183527) q[0];
x q[1];
rz(1.564304) q[2];
sx q[2];
rz(-2.3320856) q[2];
sx q[2];
rz(-0.49525515) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68542504) q[1];
sx q[1];
rz(-1.9677094) q[1];
sx q[1];
rz(-2.6422068) q[1];
rz(-1.7117483) q[3];
sx q[3];
rz(-0.8356072) q[3];
sx q[3];
rz(-0.47675374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63742796) q[2];
sx q[2];
rz(-0.97680682) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(2.6484683) q[3];
sx q[3];
rz(-0.92032856) q[3];
sx q[3];
rz(-1.5776618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.2757932) q[0];
sx q[0];
rz(-1.2340622) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(1.5419143) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(-2.7154198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0876678) q[0];
sx q[0];
rz(-3.1172981) q[0];
sx q[0];
rz(-2.3870275) q[0];
rz(-pi) q[1];
rz(2.4075899) q[2];
sx q[2];
rz(-2.2224769) q[2];
sx q[2];
rz(0.82254788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81282367) q[1];
sx q[1];
rz(-1.4534833) q[1];
sx q[1];
rz(-1.3552865) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3506463) q[3];
sx q[3];
rz(-0.90715796) q[3];
sx q[3];
rz(-2.5928465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49491945) q[2];
sx q[2];
rz(-2.9743331) q[2];
sx q[2];
rz(1.0529998) q[2];
rz(-2.7827175) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(2.0487962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36295715) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(-0.92078513) q[0];
rz(0.50865632) q[1];
sx q[1];
rz(-2.3743036) q[1];
sx q[1];
rz(2.7210534) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8201308) q[0];
sx q[0];
rz(-0.88867868) q[0];
sx q[0];
rz(0.90773037) q[0];
rz(-pi) q[1];
rz(-0.24193544) q[2];
sx q[2];
rz(-0.95429776) q[2];
sx q[2];
rz(2.7665507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9892726) q[1];
sx q[1];
rz(-2.1713683) q[1];
sx q[1];
rz(2.5737073) q[1];
rz(-2.2617662) q[3];
sx q[3];
rz(-1.4925957) q[3];
sx q[3];
rz(1.9281338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4113808) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(0.31663695) q[2];
rz(-0.5591048) q[3];
sx q[3];
rz(-0.63566256) q[3];
sx q[3];
rz(2.3234698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6560787) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(-0.47646933) q[1];
sx q[1];
rz(-2.10119) q[1];
sx q[1];
rz(1.4269921) q[1];
rz(1.477735) q[2];
sx q[2];
rz(-2.3882967) q[2];
sx q[2];
rz(-0.49327539) q[2];
rz(-1.3415402) q[3];
sx q[3];
rz(-2.7490902) q[3];
sx q[3];
rz(0.31828087) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
