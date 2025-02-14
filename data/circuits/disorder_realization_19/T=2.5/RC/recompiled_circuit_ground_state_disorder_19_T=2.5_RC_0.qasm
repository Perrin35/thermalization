OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(1.5751155) q[0];
sx q[0];
rz(10.549904) q[0];
rz(1.0220802) q[1];
sx q[1];
rz(-0.66749579) q[1];
sx q[1];
rz(-2.3281085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33685499) q[0];
sx q[0];
rz(-1.9362402) q[0];
sx q[0];
rz(0.15113896) q[0];
rz(-0.44806077) q[2];
sx q[2];
rz(-1.2431113) q[2];
sx q[2];
rz(-2.0304012) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91650016) q[1];
sx q[1];
rz(-1.4218907) q[1];
sx q[1];
rz(2.6970106) q[1];
rz(-1.7343821) q[3];
sx q[3];
rz(-1.6601052) q[3];
sx q[3];
rz(2.341604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.4514613) q[2];
sx q[2];
rz(-2.2129464) q[2];
rz(1.5993902) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(1.9699875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9376675) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(0.12705886) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(0.7712706) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.49202) q[0];
sx q[0];
rz(-2.3207773) q[0];
sx q[0];
rz(1.0787021) q[0];
x q[1];
rz(1.3992869) q[2];
sx q[2];
rz(-2.7993188) q[2];
sx q[2];
rz(-2.8245087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2576372) q[1];
sx q[1];
rz(-1.9169382) q[1];
sx q[1];
rz(2.7621072) q[1];
rz(-pi) q[2];
rz(2.1089606) q[3];
sx q[3];
rz(-1.0618883) q[3];
sx q[3];
rz(2.5050688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1473006) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(0.0083943923) q[2];
rz(-2.4781135) q[3];
sx q[3];
rz(-1.9050262) q[3];
sx q[3];
rz(-2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10107772) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(0.44152942) q[0];
rz(-0.98835522) q[1];
sx q[1];
rz(-1.1343196) q[1];
sx q[1];
rz(3.006014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5268742) q[0];
sx q[0];
rz(-1.563094) q[0];
sx q[0];
rz(3.1253424) q[0];
x q[1];
rz(-0.71696059) q[2];
sx q[2];
rz(-1.8646984) q[2];
sx q[2];
rz(0.64168054) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6511894) q[1];
sx q[1];
rz(-0.62593834) q[1];
sx q[1];
rz(-0.93749009) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3531923) q[3];
sx q[3];
rz(-0.74624589) q[3];
sx q[3];
rz(1.3775795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2077937) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(-0.03820339) q[2];
rz(2.6162052) q[3];
sx q[3];
rz(-2.5140258) q[3];
sx q[3];
rz(-2.4562522) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4811089) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(0.61087459) q[0];
rz(1.8065709) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(-0.23922051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0804028) q[0];
sx q[0];
rz(-2.1227269) q[0];
sx q[0];
rz(0.30776382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1785422) q[2];
sx q[2];
rz(-1.129727) q[2];
sx q[2];
rz(-0.46229306) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.661631) q[1];
sx q[1];
rz(-2.0212272) q[1];
sx q[1];
rz(-1.9305139) q[1];
rz(0.46142205) q[3];
sx q[3];
rz(-1.3909826) q[3];
sx q[3];
rz(1.5195886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7188344) q[2];
sx q[2];
rz(-1.6966635) q[2];
sx q[2];
rz(0.0017496721) q[2];
rz(2.8478029) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(2.8747115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9002429) q[0];
sx q[0];
rz(-1.7647864) q[0];
sx q[0];
rz(2.2985261) q[0];
rz(0.33755606) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(1.6036124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46066901) q[0];
sx q[0];
rz(-2.3516555) q[0];
sx q[0];
rz(-1.3000751) q[0];
x q[1];
rz(0.29422167) q[2];
sx q[2];
rz(-0.18202848) q[2];
sx q[2];
rz(-0.95532571) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1142769) q[1];
sx q[1];
rz(-1.3754579) q[1];
sx q[1];
rz(2.8802425) q[1];
rz(-2.9318453) q[3];
sx q[3];
rz(-1.6496067) q[3];
sx q[3];
rz(2.1926853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7097077) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(1.0464) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(-2.5939202) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043561291) q[0];
sx q[0];
rz(-1.2579608) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(-1.7421534) q[1];
sx q[1];
rz(-1.6439227) q[1];
sx q[1];
rz(1.8125777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32008313) q[0];
sx q[0];
rz(-1.3504986) q[0];
sx q[0];
rz(-3.0295323) q[0];
rz(-pi) q[1];
rz(-2.5115254) q[2];
sx q[2];
rz(-0.35379256) q[2];
sx q[2];
rz(1.9513064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4635515) q[1];
sx q[1];
rz(-1.2604509) q[1];
sx q[1];
rz(0.84872679) q[1];
x q[2];
rz(-0.51402199) q[3];
sx q[3];
rz(-1.8933305) q[3];
sx q[3];
rz(-3.1399825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38137388) q[2];
sx q[2];
rz(-1.2083283) q[2];
sx q[2];
rz(1.5927429) q[2];
rz(3.0717487) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(0.86863345) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233362) q[0];
sx q[0];
rz(-1.2160439) q[0];
sx q[0];
rz(-0.94183952) q[0];
rz(-0.16009227) q[1];
sx q[1];
rz(-1.6611049) q[1];
sx q[1];
rz(-0.19217415) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20142444) q[0];
sx q[0];
rz(-1.7946825) q[0];
sx q[0];
rz(3.1046449) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0490369) q[2];
sx q[2];
rz(-1.6055487) q[2];
sx q[2];
rz(-3.1086189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1035489) q[1];
sx q[1];
rz(-0.86665857) q[1];
sx q[1];
rz(2.0226993) q[1];
rz(-1.5302079) q[3];
sx q[3];
rz(-0.49792624) q[3];
sx q[3];
rz(2.9663309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75752246) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(0.43357098) q[2];
rz(1.5571669) q[3];
sx q[3];
rz(-1.4945364) q[3];
sx q[3];
rz(-2.7092194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(-1.3442511) q[0];
rz(1.2035707) q[1];
sx q[1];
rz(-1.1886339) q[1];
sx q[1];
rz(1.3628091) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8289611) q[0];
sx q[0];
rz(-1.5927218) q[0];
sx q[0];
rz(-1.6005105) q[0];
rz(-pi) q[1];
rz(1.7596471) q[2];
sx q[2];
rz(-0.7536234) q[2];
sx q[2];
rz(-3.0830887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92719936) q[1];
sx q[1];
rz(-1.3424557) q[1];
sx q[1];
rz(1.6990927) q[1];
rz(-pi) q[2];
rz(1.9523432) q[3];
sx q[3];
rz(-0.7233215) q[3];
sx q[3];
rz(-1.9017232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56889304) q[2];
sx q[2];
rz(-0.77304825) q[2];
sx q[2];
rz(2.1577238) q[2];
rz(-0.24108663) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-2.6702061) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(0.27467003) q[0];
rz(1.7999016) q[1];
sx q[1];
rz(-2.5119753) q[1];
sx q[1];
rz(1.0460269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5513902) q[0];
sx q[0];
rz(-2.1596163) q[0];
sx q[0];
rz(-1.9639652) q[0];
rz(-pi) q[1];
rz(-1.7685031) q[2];
sx q[2];
rz(-1.6501763) q[2];
sx q[2];
rz(-0.52620906) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1362125) q[1];
sx q[1];
rz(-2.5544205) q[1];
sx q[1];
rz(-0.11759742) q[1];
rz(-2.0738129) q[3];
sx q[3];
rz(-1.6075589) q[3];
sx q[3];
rz(-3.0779148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(2.7421303) q[2];
rz(-0.56600371) q[3];
sx q[3];
rz(-0.96650201) q[3];
sx q[3];
rz(-0.81542265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28867662) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(-0.5823108) q[0];
rz(-2.9583926) q[1];
sx q[1];
rz(-2.2661426) q[1];
sx q[1];
rz(-0.80642548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69003478) q[0];
sx q[0];
rz(-1.6584883) q[0];
sx q[0];
rz(1.466485) q[0];
rz(1.3034794) q[2];
sx q[2];
rz(-1.2462933) q[2];
sx q[2];
rz(-2.7410206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.90102531) q[1];
sx q[1];
rz(-1.7738924) q[1];
sx q[1];
rz(0.74764772) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6819477) q[3];
sx q[3];
rz(-1.2254997) q[3];
sx q[3];
rz(0.87792859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9091984) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(-2.3363028) q[2];
rz(-2.9946839) q[3];
sx q[3];
rz(-0.78607905) q[3];
sx q[3];
rz(-2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2522226) q[0];
sx q[0];
rz(-1.7563553) q[0];
sx q[0];
rz(1.8808543) q[0];
rz(-1.5994785) q[1];
sx q[1];
rz(-1.5442994) q[1];
sx q[1];
rz(-1.50179) q[1];
rz(3.0270544) q[2];
sx q[2];
rz(-0.85083148) q[2];
sx q[2];
rz(1.1640805) q[2];
rz(1.2272502) q[3];
sx q[3];
rz(-0.95442964) q[3];
sx q[3];
rz(-1.7908532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
