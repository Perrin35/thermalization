OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6516946) q[0];
sx q[0];
rz(5.4013847) q[0];
sx q[0];
rz(9.2800979) q[0];
rz(-2.7235003) q[1];
sx q[1];
rz(-3.0378208) q[1];
sx q[1];
rz(1.6296847) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1967752) q[0];
sx q[0];
rz(-2.631105) q[0];
sx q[0];
rz(1.8148242) q[0];
x q[1];
rz(-1.6250526) q[2];
sx q[2];
rz(-2.3035604) q[2];
sx q[2];
rz(-0.65544477) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.50070333) q[1];
sx q[1];
rz(-1.902406) q[1];
sx q[1];
rz(-1.9808116) q[1];
rz(0.6717967) q[3];
sx q[3];
rz(-0.62272391) q[3];
sx q[3];
rz(-1.1237564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8251553) q[2];
sx q[2];
rz(-1.8708159) q[2];
sx q[2];
rz(-2.4862508) q[2];
rz(-1.7573028) q[3];
sx q[3];
rz(-0.96450788) q[3];
sx q[3];
rz(2.3206319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69219387) q[0];
sx q[0];
rz(-1.7658424) q[0];
sx q[0];
rz(2.6334515) q[0];
rz(1.7610158) q[1];
sx q[1];
rz(-2.5602129) q[1];
sx q[1];
rz(1.2095721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66099453) q[0];
sx q[0];
rz(-0.61836243) q[0];
sx q[0];
rz(-1.2510651) q[0];
rz(-pi) q[1];
rz(-1.5530994) q[2];
sx q[2];
rz(-1.624384) q[2];
sx q[2];
rz(-2.686116) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0703549) q[1];
sx q[1];
rz(-2.6609328) q[1];
sx q[1];
rz(2.8621469) q[1];
rz(-pi) q[2];
rz(0.38755667) q[3];
sx q[3];
rz(-1.0021375) q[3];
sx q[3];
rz(-2.89944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2035344) q[2];
sx q[2];
rz(-1.396023) q[2];
sx q[2];
rz(0.61689845) q[2];
rz(1.6054224) q[3];
sx q[3];
rz(-1.0441531) q[3];
sx q[3];
rz(2.6507586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9817552) q[0];
sx q[0];
rz(-1.6438537) q[0];
sx q[0];
rz(-0.72506654) q[0];
rz(-1.2695351) q[1];
sx q[1];
rz(-2.4639362) q[1];
sx q[1];
rz(-2.1824172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77391439) q[0];
sx q[0];
rz(-1.1836021) q[0];
sx q[0];
rz(-1.0605124) q[0];
x q[1];
rz(-1.959402) q[2];
sx q[2];
rz(-1.1827785) q[2];
sx q[2];
rz(1.6418552) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.945248) q[1];
sx q[1];
rz(-1.3875393) q[1];
sx q[1];
rz(3.0250579) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3884344) q[3];
sx q[3];
rz(-1.404247) q[3];
sx q[3];
rz(-0.20660755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3265257) q[2];
sx q[2];
rz(-1.8748137) q[2];
sx q[2];
rz(0.25700021) q[2];
rz(-2.0604996) q[3];
sx q[3];
rz(-2.6205781) q[3];
sx q[3];
rz(1.5787554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0839888) q[0];
sx q[0];
rz(-2.8631518) q[0];
sx q[0];
rz(1.1778911) q[0];
rz(-0.15011694) q[1];
sx q[1];
rz(-0.53214407) q[1];
sx q[1];
rz(-1.5163126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0310136) q[0];
sx q[0];
rz(-0.4415126) q[0];
sx q[0];
rz(-0.86236055) q[0];
rz(-0.37090918) q[2];
sx q[2];
rz(-2.264655) q[2];
sx q[2];
rz(2.3313076) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2840239) q[1];
sx q[1];
rz(-2.800436) q[1];
sx q[1];
rz(-0.30156231) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4978906) q[3];
sx q[3];
rz(-1.4182404) q[3];
sx q[3];
rz(1.2920472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2991422) q[2];
sx q[2];
rz(-2.5717042) q[2];
sx q[2];
rz(2.2184856) q[2];
rz(-0.9592157) q[3];
sx q[3];
rz(-1.0126637) q[3];
sx q[3];
rz(-2.925351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32090309) q[0];
sx q[0];
rz(-0.90383363) q[0];
sx q[0];
rz(0.88291105) q[0];
rz(-1.8461022) q[1];
sx q[1];
rz(-0.77406445) q[1];
sx q[1];
rz(-0.49497089) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8710324) q[0];
sx q[0];
rz(-0.43631662) q[0];
sx q[0];
rz(-2.3113234) q[0];
rz(0.70303838) q[2];
sx q[2];
rz(-1.993317) q[2];
sx q[2];
rz(-0.51681821) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7385612) q[1];
sx q[1];
rz(-0.90442362) q[1];
sx q[1];
rz(-1.6825466) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0932156) q[3];
sx q[3];
rz(-1.256878) q[3];
sx q[3];
rz(0.81410223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8098658) q[2];
sx q[2];
rz(-0.11955424) q[2];
sx q[2];
rz(-1.8734107) q[2];
rz(0.082503334) q[3];
sx q[3];
rz(-1.5039597) q[3];
sx q[3];
rz(-2.4845607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4853915) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(-1.50151) q[0];
rz(2.0791176) q[1];
sx q[1];
rz(-1.1524009) q[1];
sx q[1];
rz(1.9662439) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1262681) q[0];
sx q[0];
rz(-1.7909174) q[0];
sx q[0];
rz(-0.31173978) q[0];
x q[1];
rz(-2.7926867) q[2];
sx q[2];
rz(-2.4173173) q[2];
sx q[2];
rz(-1.6394212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8379587) q[1];
sx q[1];
rz(-0.278618) q[1];
sx q[1];
rz(-2.7024804) q[1];
rz(-pi) q[2];
rz(2.8296521) q[3];
sx q[3];
rz(-2.3565203) q[3];
sx q[3];
rz(-1.5842445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33092734) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(-2.1539099) q[2];
rz(0.87743131) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(0.73006829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9577318) q[0];
sx q[0];
rz(-1.6022302) q[0];
sx q[0];
rz(3.0027332) q[0];
rz(1.214437) q[1];
sx q[1];
rz(-1.6338394) q[1];
sx q[1];
rz(0.81659281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3735376) q[0];
sx q[0];
rz(-1.5898503) q[0];
sx q[0];
rz(1.7956177) q[0];
rz(-0.26890484) q[2];
sx q[2];
rz(-2.3550526) q[2];
sx q[2];
rz(-1.5603017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5948171) q[1];
sx q[1];
rz(-0.71974126) q[1];
sx q[1];
rz(-1.5139346) q[1];
rz(-pi) q[2];
rz(-0.66612996) q[3];
sx q[3];
rz(-2.1762848) q[3];
sx q[3];
rz(-0.44119409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3198118) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(-0.4846586) q[2];
rz(-0.36568493) q[3];
sx q[3];
rz(-1.7771143) q[3];
sx q[3];
rz(1.9827838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1300238) q[0];
sx q[0];
rz(-0.31112177) q[0];
sx q[0];
rz(-0.097231641) q[0];
rz(-0.10487996) q[1];
sx q[1];
rz(-0.77969867) q[1];
sx q[1];
rz(-1.8269151) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0630232) q[0];
sx q[0];
rz(-0.0067575909) q[0];
sx q[0];
rz(-1.2122173) q[0];
rz(-pi) q[1];
rz(-1.329256) q[2];
sx q[2];
rz(-1.7881696) q[2];
sx q[2];
rz(-1.5474943) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2995616) q[1];
sx q[1];
rz(-2.6202469) q[1];
sx q[1];
rz(3.035665) q[1];
rz(-2.1076074) q[3];
sx q[3];
rz(-1.4958463) q[3];
sx q[3];
rz(-1.0190462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17681992) q[2];
sx q[2];
rz(-1.7440045) q[2];
sx q[2];
rz(2.9883265) q[2];
rz(-1.2166474) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(0.6161859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0230873) q[0];
sx q[0];
rz(-3.1361134) q[0];
sx q[0];
rz(-1.6368921) q[0];
rz(2.3432689) q[1];
sx q[1];
rz(-1.5232122) q[1];
sx q[1];
rz(-0.33531478) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9630747) q[0];
sx q[0];
rz(-1.3097242) q[0];
sx q[0];
rz(0.065351323) q[0];
rz(-3.1153572) q[2];
sx q[2];
rz(-2.3910284) q[2];
sx q[2];
rz(-0.98294965) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7882725) q[1];
sx q[1];
rz(-1.2306552) q[1];
sx q[1];
rz(1.2659094) q[1];
rz(-2.9842989) q[3];
sx q[3];
rz(-0.99911896) q[3];
sx q[3];
rz(-1.8974822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0140784) q[2];
sx q[2];
rz(-1.9731584) q[2];
sx q[2];
rz(-2.7039995) q[2];
rz(-2.0994999) q[3];
sx q[3];
rz(-1.7362678) q[3];
sx q[3];
rz(-2.5991345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1522778) q[0];
sx q[0];
rz(-0.88459009) q[0];
sx q[0];
rz(2.6919795) q[0];
rz(0.97688976) q[1];
sx q[1];
rz(-1.5039109) q[1];
sx q[1];
rz(2.6403715) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2758482) q[0];
sx q[0];
rz(-1.8536708) q[0];
sx q[0];
rz(-0.19794835) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0222763) q[2];
sx q[2];
rz(-0.46487936) q[2];
sx q[2];
rz(1.3550188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38793716) q[1];
sx q[1];
rz(-0.5273312) q[1];
sx q[1];
rz(2.0547632) q[1];
rz(-0.95285033) q[3];
sx q[3];
rz(-2.5898159) q[3];
sx q[3];
rz(-2.1976041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0235128) q[2];
sx q[2];
rz(-0.89236516) q[2];
sx q[2];
rz(-0.21772131) q[2];
rz(1.2348385) q[3];
sx q[3];
rz(-2.4775938) q[3];
sx q[3];
rz(2.2209404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6476743) q[0];
sx q[0];
rz(-1.7128581) q[0];
sx q[0];
rz(0.38059522) q[0];
rz(0.71161288) q[1];
sx q[1];
rz(-1.8743534) q[1];
sx q[1];
rz(1.4030917) q[1];
rz(-3.0416476) q[2];
sx q[2];
rz(-2.6700085) q[2];
sx q[2];
rz(2.811583) q[2];
rz(2.6307132) q[3];
sx q[3];
rz(-2.6838959) q[3];
sx q[3];
rz(0.30797227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
