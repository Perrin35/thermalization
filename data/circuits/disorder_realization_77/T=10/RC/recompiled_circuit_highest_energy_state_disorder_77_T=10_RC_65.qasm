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
rz(-0.51802975) q[0];
sx q[0];
rz(-2.0002444) q[0];
sx q[0];
rz(-0.040123392) q[0];
rz(-2.8934381) q[1];
sx q[1];
rz(-2.0791972) q[1];
sx q[1];
rz(0.084029347) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2618128) q[0];
sx q[0];
rz(-0.09238681) q[0];
sx q[0];
rz(-2.4506344) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3925926) q[2];
sx q[2];
rz(-1.1786818) q[2];
sx q[2];
rz(-0.41022476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58761105) q[1];
sx q[1];
rz(-1.1003255) q[1];
sx q[1];
rz(1.379983) q[1];
rz(-1.1575451) q[3];
sx q[3];
rz(-1.0344369) q[3];
sx q[3];
rz(-2.9265917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4485126) q[2];
sx q[2];
rz(-1.7757519) q[2];
sx q[2];
rz(-0.25126323) q[2];
rz(-1.9692028) q[3];
sx q[3];
rz(-1.1029714) q[3];
sx q[3];
rz(0.049886543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6867111) q[0];
sx q[0];
rz(-2.659681) q[0];
sx q[0];
rz(1.7676109) q[0];
rz(0.33085597) q[1];
sx q[1];
rz(-1.4127981) q[1];
sx q[1];
rz(0.27423283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62456709) q[0];
sx q[0];
rz(-0.32450482) q[0];
sx q[0];
rz(-0.22402482) q[0];
rz(-pi) q[1];
rz(-2.6083228) q[2];
sx q[2];
rz(-1.5646213) q[2];
sx q[2];
rz(1.285163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9668558) q[1];
sx q[1];
rz(-1.2464332) q[1];
sx q[1];
rz(-1.3423389) q[1];
rz(-1.7555439) q[3];
sx q[3];
rz(-2.7165301) q[3];
sx q[3];
rz(-2.5109072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26068822) q[2];
sx q[2];
rz(-1.8131249) q[2];
sx q[2];
rz(-2.0534959) q[2];
rz(2.0982096) q[3];
sx q[3];
rz(-2.9731396) q[3];
sx q[3];
rz(-2.8823631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.008721) q[0];
sx q[0];
rz(-2.6368124) q[0];
sx q[0];
rz(1.1016499) q[0];
rz(2.3654826) q[1];
sx q[1];
rz(-1.2932237) q[1];
sx q[1];
rz(-0.64220846) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4896079) q[0];
sx q[0];
rz(-2.5074158) q[0];
sx q[0];
rz(-0.80192566) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71977653) q[2];
sx q[2];
rz(-1.7584561) q[2];
sx q[2];
rz(-1.3087147) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5635665) q[1];
sx q[1];
rz(-1.9959004) q[1];
sx q[1];
rz(-0.11386392) q[1];
rz(-2.4829743) q[3];
sx q[3];
rz(-2.082654) q[3];
sx q[3];
rz(-3.1114674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3099826) q[2];
sx q[2];
rz(-1.1534561) q[2];
sx q[2];
rz(2.2115121) q[2];
rz(-2.3275404) q[3];
sx q[3];
rz(-1.1907153) q[3];
sx q[3];
rz(-1.2306151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.880068) q[0];
sx q[0];
rz(-1.1898758) q[0];
sx q[0];
rz(2.6370908) q[0];
rz(-2.2241459) q[1];
sx q[1];
rz(-2.4023299) q[1];
sx q[1];
rz(-2.9333072) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2613647) q[0];
sx q[0];
rz(-2.8552736) q[0];
sx q[0];
rz(2.5891735) q[0];
rz(0.31625749) q[2];
sx q[2];
rz(-2.3964747) q[2];
sx q[2];
rz(0.032501246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5985344) q[1];
sx q[1];
rz(-1.0987135) q[1];
sx q[1];
rz(0.92453875) q[1];
rz(-pi) q[2];
rz(-2.743163) q[3];
sx q[3];
rz(-2.5442985) q[3];
sx q[3];
rz(0.78890991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31671277) q[2];
sx q[2];
rz(-1.0944288) q[2];
sx q[2];
rz(1.1701976) q[2];
rz(-0.96585387) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(2.9962208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.589094) q[0];
sx q[0];
rz(-2.4674802) q[0];
sx q[0];
rz(0.66459769) q[0];
rz(0.26847863) q[1];
sx q[1];
rz(-1.6842027) q[1];
sx q[1];
rz(2.1427515) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12832175) q[0];
sx q[0];
rz(-1.7230464) q[0];
sx q[0];
rz(0.98441846) q[0];
x q[1];
rz(1.8605609) q[2];
sx q[2];
rz(-2.701607) q[2];
sx q[2];
rz(2.5948465) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4722574) q[1];
sx q[1];
rz(-0.20755033) q[1];
sx q[1];
rz(2.186354) q[1];
rz(-pi) q[2];
rz(-0.020962997) q[3];
sx q[3];
rz(-0.87910324) q[3];
sx q[3];
rz(-1.5251336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9773679) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(-0.11717907) q[2];
rz(3.0818648) q[3];
sx q[3];
rz(-1.9220587) q[3];
sx q[3];
rz(0.81502771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77434671) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(1.8814948) q[0];
rz(-3.0267808) q[1];
sx q[1];
rz(-1.7962619) q[1];
sx q[1];
rz(3.0462435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0336825) q[0];
sx q[0];
rz(-2.9386407) q[0];
sx q[0];
rz(0.072921948) q[0];
rz(-pi) q[1];
rz(-0.16654715) q[2];
sx q[2];
rz(-1.0726561) q[2];
sx q[2];
rz(2.0280251) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29987511) q[1];
sx q[1];
rz(-2.5073194) q[1];
sx q[1];
rz(0.81171616) q[1];
x q[2];
rz(3.1308858) q[3];
sx q[3];
rz(-2.0981776) q[3];
sx q[3];
rz(0.37400613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.054333869) q[2];
sx q[2];
rz(-2.2237873) q[2];
sx q[2];
rz(2.4799662) q[2];
rz(-2.3719487) q[3];
sx q[3];
rz(-1.6911643) q[3];
sx q[3];
rz(2.2746287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041644) q[0];
sx q[0];
rz(-0.34808174) q[0];
sx q[0];
rz(0.75781703) q[0];
rz(0.025253145) q[1];
sx q[1];
rz(-0.39552894) q[1];
sx q[1];
rz(1.1753561) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379743) q[0];
sx q[0];
rz(-1.0991862) q[0];
sx q[0];
rz(-1.1704117) q[0];
x q[1];
rz(-0.98837672) q[2];
sx q[2];
rz(-2.3727131) q[2];
sx q[2];
rz(-2.0748367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7832406) q[1];
sx q[1];
rz(-2.2951627) q[1];
sx q[1];
rz(1.0388166) q[1];
x q[2];
rz(1.5478915) q[3];
sx q[3];
rz(-1.1419019) q[3];
sx q[3];
rz(-1.0064841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5777099) q[2];
sx q[2];
rz(-1.5241728) q[2];
sx q[2];
rz(1.9263402) q[2];
rz(-0.75872129) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(-1.5846213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6587081) q[0];
sx q[0];
rz(-1.8526798) q[0];
sx q[0];
rz(-0.38920745) q[0];
rz(3.0534741) q[1];
sx q[1];
rz(-1.7033109) q[1];
sx q[1];
rz(-1.6099991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22613285) q[0];
sx q[0];
rz(-2.8185039) q[0];
sx q[0];
rz(0.74193546) q[0];
rz(-0.66122239) q[2];
sx q[2];
rz(-1.9384137) q[2];
sx q[2];
rz(-1.6501381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1080146) q[1];
sx q[1];
rz(-0.44541767) q[1];
sx q[1];
rz(-0.12081318) q[1];
rz(-2.0069028) q[3];
sx q[3];
rz(-2.3824771) q[3];
sx q[3];
rz(1.5576943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4072998) q[2];
sx q[2];
rz(-0.6137085) q[2];
sx q[2];
rz(-1.8180004) q[2];
rz(1.7981516) q[3];
sx q[3];
rz(-0.7898134) q[3];
sx q[3];
rz(2.8048803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6650498) q[0];
sx q[0];
rz(-0.90710586) q[0];
sx q[0];
rz(-0.46515775) q[0];
rz(-3.0910885) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(-2.6506298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28117958) q[0];
sx q[0];
rz(-1.2523012) q[0];
sx q[0];
rz(2.7943576) q[0];
rz(-pi) q[1];
rz(-0.45456205) q[2];
sx q[2];
rz(-2.801932) q[2];
sx q[2];
rz(-2.8220334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9742187) q[1];
sx q[1];
rz(-0.95229665) q[1];
sx q[1];
rz(2.2688686) q[1];
x q[2];
rz(-2.773953) q[3];
sx q[3];
rz(-0.30045569) q[3];
sx q[3];
rz(-0.57747546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65840536) q[2];
sx q[2];
rz(-0.5961954) q[2];
sx q[2];
rz(2.2779706) q[2];
rz(-2.9577799) q[3];
sx q[3];
rz(-2.176216) q[3];
sx q[3];
rz(0.98627728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.42613906) q[0];
sx q[0];
rz(-1.5197536) q[0];
sx q[0];
rz(2.5157628) q[0];
rz(-0.65678701) q[1];
sx q[1];
rz(-0.96581179) q[1];
sx q[1];
rz(1.2084557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0330808) q[0];
sx q[0];
rz(-1.4377366) q[0];
sx q[0];
rz(3.0390374) q[0];
rz(1.8911458) q[2];
sx q[2];
rz(-0.9365754) q[2];
sx q[2];
rz(-2.5958217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.085026) q[1];
sx q[1];
rz(-1.4843656) q[1];
sx q[1];
rz(1.568145) q[1];
rz(1.7478862) q[3];
sx q[3];
rz(-2.13518) q[3];
sx q[3];
rz(-2.0616814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5586231) q[2];
sx q[2];
rz(-1.308459) q[2];
sx q[2];
rz(2.6226131) q[2];
rz(2.6988622) q[3];
sx q[3];
rz(-1.818592) q[3];
sx q[3];
rz(-1.3435266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2124355) q[0];
sx q[0];
rz(-0.99016187) q[0];
sx q[0];
rz(-1.1708175) q[0];
rz(-0.15405542) q[1];
sx q[1];
rz(-0.93692056) q[1];
sx q[1];
rz(-0.9137203) q[1];
rz(-1.6120934) q[2];
sx q[2];
rz(-2.0870024) q[2];
sx q[2];
rz(0.11556297) q[2];
rz(-0.89755015) q[3];
sx q[3];
rz(-1.8445476) q[3];
sx q[3];
rz(-2.6060819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
