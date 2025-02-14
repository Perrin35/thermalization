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
rz(-2.2588377) q[0];
sx q[0];
rz(-0.14482276) q[0];
sx q[0];
rz(2.3902399) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(5.3137988) q[1];
sx q[1];
rz(9.245524) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0517901) q[0];
sx q[0];
rz(-1.8500916) q[0];
sx q[0];
rz(-1.0599049) q[0];
rz(-pi) q[1];
rz(-1.5380074) q[2];
sx q[2];
rz(-2.2934487) q[2];
sx q[2];
rz(3.0422831) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1826659) q[1];
sx q[1];
rz(-1.2616871) q[1];
sx q[1];
rz(1.6462417) q[1];
rz(-pi) q[2];
rz(-0.93710812) q[3];
sx q[3];
rz(-1.6185915) q[3];
sx q[3];
rz(0.60768793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2396607) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(0.61061668) q[2];
rz(-0.88879746) q[3];
sx q[3];
rz(-2.461268) q[3];
sx q[3];
rz(-2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4438542) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(1.8119716) q[0];
rz(-1.4854206) q[1];
sx q[1];
rz(-1.679436) q[1];
sx q[1];
rz(-2.2775473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5764424) q[0];
sx q[0];
rz(-0.2731384) q[0];
sx q[0];
rz(-1.1821724) q[0];
rz(2.0338221) q[2];
sx q[2];
rz(-1.8512176) q[2];
sx q[2];
rz(-0.011649557) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4046574) q[1];
sx q[1];
rz(-2.9184249) q[1];
sx q[1];
rz(-2.2523802) q[1];
x q[2];
rz(-2.1883499) q[3];
sx q[3];
rz(-2.2526178) q[3];
sx q[3];
rz(-0.98495959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0835421) q[2];
sx q[2];
rz(-0.6821878) q[2];
sx q[2];
rz(0.19134276) q[2];
rz(-2.8355016) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3092344) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(-0.424463) q[0];
rz(-1.8008495) q[1];
sx q[1];
rz(-0.91573358) q[1];
sx q[1];
rz(-0.65139604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0661591) q[0];
sx q[0];
rz(-1.4071583) q[0];
sx q[0];
rz(3.1343824) q[0];
rz(-pi) q[1];
rz(-0.38369757) q[2];
sx q[2];
rz(-1.7406775) q[2];
sx q[2];
rz(3.1184514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9561852) q[1];
sx q[1];
rz(-0.53403097) q[1];
sx q[1];
rz(2.7772481) q[1];
x q[2];
rz(-0.15768361) q[3];
sx q[3];
rz(-1.3467977) q[3];
sx q[3];
rz(1.7832613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1059619) q[2];
sx q[2];
rz(-1.1392081) q[2];
sx q[2];
rz(-2.4448709) q[2];
rz(2.4524955) q[3];
sx q[3];
rz(-1.9870116) q[3];
sx q[3];
rz(2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43831393) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(-1.0003723) q[0];
rz(-1.6814303) q[1];
sx q[1];
rz(-1.5337475) q[1];
sx q[1];
rz(-1.9283074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.723169) q[0];
sx q[0];
rz(-0.30659404) q[0];
sx q[0];
rz(-1.1585537) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4126265) q[2];
sx q[2];
rz(-1.2104958) q[2];
sx q[2];
rz(-0.0055706669) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91649969) q[1];
sx q[1];
rz(-1.8986393) q[1];
sx q[1];
rz(1.4429379) q[1];
rz(-0.090661006) q[3];
sx q[3];
rz(-2.1075776) q[3];
sx q[3];
rz(1.2556374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0393684) q[2];
sx q[2];
rz(-2.3185456) q[2];
sx q[2];
rz(-3.0774806) q[2];
rz(-0.72470775) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0475912) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(0.79175788) q[0];
rz(-2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2091141) q[0];
sx q[0];
rz(-0.9830342) q[0];
sx q[0];
rz(0.57508075) q[0];
rz(0.58081268) q[2];
sx q[2];
rz(-1.6585095) q[2];
sx q[2];
rz(0.49077362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2446652) q[1];
sx q[1];
rz(-1.7411971) q[1];
sx q[1];
rz(-0.076431304) q[1];
rz(0.20466699) q[3];
sx q[3];
rz(-1.0042666) q[3];
sx q[3];
rz(0.03290225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84805924) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(-1.1406356) q[2];
rz(3.0349777) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(2.831736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3265729) q[0];
sx q[0];
rz(-2.6928379) q[0];
sx q[0];
rz(3.1383681) q[0];
rz(0.047317304) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(1.7914194) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2285706) q[0];
sx q[0];
rz(-1.6294894) q[0];
sx q[0];
rz(1.8236158) q[0];
x q[1];
rz(-2.2590738) q[2];
sx q[2];
rz(-0.79791683) q[2];
sx q[2];
rz(1.5882815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1941095) q[1];
sx q[1];
rz(-1.0110185) q[1];
sx q[1];
rz(-2.4653788) q[1];
x q[2];
rz(0.030478625) q[3];
sx q[3];
rz(-1.7913831) q[3];
sx q[3];
rz(0.36530803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1515767) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(3.0604176) q[2];
rz(-2.2422527) q[3];
sx q[3];
rz(-2.3266561) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154295) q[0];
sx q[0];
rz(-1.8303215) q[0];
sx q[0];
rz(2.6370866) q[0];
rz(-0.99507487) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(0.02034932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58344719) q[0];
sx q[0];
rz(-0.37675315) q[0];
sx q[0];
rz(1.3135629) q[0];
x q[1];
rz(3.0145524) q[2];
sx q[2];
rz(-1.7774701) q[2];
sx q[2];
rz(-0.57422306) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87096244) q[1];
sx q[1];
rz(-0.32901627) q[1];
sx q[1];
rz(0.94493072) q[1];
x q[2];
rz(1.0579823) q[3];
sx q[3];
rz(-0.098976299) q[3];
sx q[3];
rz(1.3266226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67123479) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(-1.7669558) q[2];
rz(-1.6124407) q[3];
sx q[3];
rz(-2.0917454) q[3];
sx q[3];
rz(3.0488739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2328211) q[0];
sx q[0];
rz(-3.0290373) q[0];
sx q[0];
rz(-2.0594647) q[0];
rz(0.96083653) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(0.94295162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0531702) q[0];
sx q[0];
rz(-2.5665847) q[0];
sx q[0];
rz(-1.0002329) q[0];
x q[1];
rz(-0.93859886) q[2];
sx q[2];
rz(-1.7603121) q[2];
sx q[2];
rz(1.0908443) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8170191) q[1];
sx q[1];
rz(-0.96265154) q[1];
sx q[1];
rz(0.84546802) q[1];
x q[2];
rz(1.2528768) q[3];
sx q[3];
rz(-2.8189481) q[3];
sx q[3];
rz(1.2979729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97041398) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(1.1478434) q[2];
rz(2.7796699) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(-1.743478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.73087937) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(-3.0391589) q[0];
rz(-2.5275285) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(-2.3238497) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2376643) q[0];
sx q[0];
rz(-0.92887628) q[0];
sx q[0];
rz(-0.68591811) q[0];
x q[1];
rz(2.2417415) q[2];
sx q[2];
rz(-1.6616115) q[2];
sx q[2];
rz(-0.95041529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2008668) q[1];
sx q[1];
rz(-1.5897017) q[1];
sx q[1];
rz(-2.3000642) q[1];
x q[2];
rz(-2.1675112) q[3];
sx q[3];
rz(-0.72457216) q[3];
sx q[3];
rz(-0.71648471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7927336) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(2.8279772) q[2];
rz(2.6775728) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(-0.919842) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87026507) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(-0.23649293) q[0];
rz(2.927921) q[1];
sx q[1];
rz(-2.2387319) q[1];
sx q[1];
rz(-0.22458354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29447039) q[0];
sx q[0];
rz(-1.0884388) q[0];
sx q[0];
rz(-0.76089528) q[0];
rz(0.5245536) q[2];
sx q[2];
rz(-0.65215014) q[2];
sx q[2];
rz(2.1422276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.869369) q[1];
sx q[1];
rz(-1.0450796) q[1];
sx q[1];
rz(1.7030667) q[1];
rz(-pi) q[2];
rz(0.41992374) q[3];
sx q[3];
rz(-2.2444199) q[3];
sx q[3];
rz(-0.47823423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.97682041) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(0.15178794) q[2];
rz(3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(3.013124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0781773) q[0];
sx q[0];
rz(-1.6178394) q[0];
sx q[0];
rz(2.7375426) q[0];
rz(-2.0582485) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(0.2486817) q[2];
sx q[2];
rz(-2.9338825) q[2];
sx q[2];
rz(2.6354229) q[2];
rz(-1.6733698) q[3];
sx q[3];
rz(-1.8986618) q[3];
sx q[3];
rz(-0.31019216) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
