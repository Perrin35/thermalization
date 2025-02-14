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
rz(-3.1168923) q[0];
sx q[0];
rz(-1.8827117) q[0];
sx q[0];
rz(-1.2727241) q[0];
rz(1.3920353) q[1];
sx q[1];
rz(-1.4104383) q[1];
sx q[1];
rz(-0.8937723) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05554219) q[0];
sx q[0];
rz(-2.1889733) q[0];
sx q[0];
rz(2.2937141) q[0];
x q[1];
rz(-0.29869334) q[2];
sx q[2];
rz(-1.8161886) q[2];
sx q[2];
rz(1.1358062) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3113678) q[1];
sx q[1];
rz(-2.0067978) q[1];
sx q[1];
rz(0.59200759) q[1];
x q[2];
rz(3.0953206) q[3];
sx q[3];
rz(-0.36847575) q[3];
sx q[3];
rz(-0.61532332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7749403) q[2];
sx q[2];
rz(-1.4619991) q[2];
sx q[2];
rz(-0.53202638) q[2];
rz(-2.4074647) q[3];
sx q[3];
rz(-2.8668154) q[3];
sx q[3];
rz(-1.1067357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72320139) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(-0.080408737) q[0];
rz(-0.53781992) q[1];
sx q[1];
rz(-1.0686914) q[1];
sx q[1];
rz(0.72371662) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2576669) q[0];
sx q[0];
rz(-2.2307255) q[0];
sx q[0];
rz(-2.1377449) q[0];
x q[1];
rz(2.2494456) q[2];
sx q[2];
rz(-0.79024678) q[2];
sx q[2];
rz(0.21871601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.58077512) q[1];
sx q[1];
rz(-2.5530949) q[1];
sx q[1];
rz(1.5735516) q[1];
rz(-pi) q[2];
rz(-0.51809727) q[3];
sx q[3];
rz(-0.88740098) q[3];
sx q[3];
rz(0.70113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2878652) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(2.2354324) q[2];
rz(-0.3624889) q[3];
sx q[3];
rz(-1.5443085) q[3];
sx q[3];
rz(0.046886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2694117) q[0];
sx q[0];
rz(-1.8159287) q[0];
sx q[0];
rz(1.0915225) q[0];
rz(-3.0190234) q[1];
sx q[1];
rz(-1.4361607) q[1];
sx q[1];
rz(-1.3465808) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4849834) q[0];
sx q[0];
rz(-0.35728982) q[0];
sx q[0];
rz(1.6553418) q[0];
rz(-pi) q[1];
rz(2.7287219) q[2];
sx q[2];
rz(-0.70253583) q[2];
sx q[2];
rz(-1.2324049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8217762) q[1];
sx q[1];
rz(-0.21142658) q[1];
sx q[1];
rz(-1.3216622) q[1];
x q[2];
rz(-2.9503959) q[3];
sx q[3];
rz(-0.53732291) q[3];
sx q[3];
rz(2.4464698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71318212) q[2];
sx q[2];
rz(-2.0117663) q[2];
sx q[2];
rz(2.909519) q[2];
rz(-1.9706767) q[3];
sx q[3];
rz(-0.62058312) q[3];
sx q[3];
rz(2.8569729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93037) q[0];
sx q[0];
rz(-0.79293293) q[0];
sx q[0];
rz(1.5027745) q[0];
rz(-0.59208313) q[1];
sx q[1];
rz(-1.1555669) q[1];
sx q[1];
rz(0.88964644) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9735721) q[0];
sx q[0];
rz(-2.4102927) q[0];
sx q[0];
rz(0.84748737) q[0];
rz(2.0621262) q[2];
sx q[2];
rz(-1.7867807) q[2];
sx q[2];
rz(0.94397533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6742176) q[1];
sx q[1];
rz(-1.8867889) q[1];
sx q[1];
rz(0.18958474) q[1];
x q[2];
rz(-0.82075228) q[3];
sx q[3];
rz(-2.6057557) q[3];
sx q[3];
rz(-1.7465296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4952116) q[2];
sx q[2];
rz(-1.6782328) q[2];
sx q[2];
rz(-0.67592534) q[2];
rz(1.4365139) q[3];
sx q[3];
rz(-2.8016475) q[3];
sx q[3];
rz(-0.16726141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2109569) q[0];
sx q[0];
rz(-0.3322596) q[0];
sx q[0];
rz(-0.19530547) q[0];
rz(-0.12061128) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(2.9511071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2094177) q[0];
sx q[0];
rz(-0.10917347) q[0];
sx q[0];
rz(-2.5072844) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.484507) q[2];
sx q[2];
rz(-0.89601529) q[2];
sx q[2];
rz(-0.23160411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6836759) q[1];
sx q[1];
rz(-1.4701587) q[1];
sx q[1];
rz(1.6229873) q[1];
rz(-pi) q[2];
rz(2.5267532) q[3];
sx q[3];
rz(-1.9214464) q[3];
sx q[3];
rz(-2.9235554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5087937) q[2];
sx q[2];
rz(-1.4845279) q[2];
sx q[2];
rz(1.9602027) q[2];
rz(1.3440291) q[3];
sx q[3];
rz(-0.7414147) q[3];
sx q[3];
rz(-0.77176362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33264273) q[0];
sx q[0];
rz(-1.3430261) q[0];
sx q[0];
rz(2.0020265) q[0];
rz(2.471916) q[1];
sx q[1];
rz(-2.2256336) q[1];
sx q[1];
rz(1.12961) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07189508) q[0];
sx q[0];
rz(-0.41145936) q[0];
sx q[0];
rz(-0.81019391) q[0];
x q[1];
rz(1.8367581) q[2];
sx q[2];
rz(-0.87087357) q[2];
sx q[2];
rz(-2.093932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.672513) q[1];
sx q[1];
rz(-1.3762849) q[1];
sx q[1];
rz(-2.7698293) q[1];
rz(-pi) q[2];
rz(-2.9608742) q[3];
sx q[3];
rz(-1.7263432) q[3];
sx q[3];
rz(-0.29744086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42664042) q[2];
sx q[2];
rz(-1.4098097) q[2];
sx q[2];
rz(-0.40531522) q[2];
rz(-0.82695588) q[3];
sx q[3];
rz(-0.6260286) q[3];
sx q[3];
rz(0.85957447) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3655132) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(-0.93589163) q[0];
rz(2.2731958) q[1];
sx q[1];
rz(-2.6135542) q[1];
sx q[1];
rz(1.7220928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7129684) q[0];
sx q[0];
rz(-2.0789062) q[0];
sx q[0];
rz(1.1626121) q[0];
x q[1];
rz(1.1677891) q[2];
sx q[2];
rz(-0.54318494) q[2];
sx q[2];
rz(-2.3576971) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4861814) q[1];
sx q[1];
rz(-2.2194194) q[1];
sx q[1];
rz(0.59452615) q[1];
x q[2];
rz(1.5711745) q[3];
sx q[3];
rz(-1.5690363) q[3];
sx q[3];
rz(1.5906217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3203656) q[2];
sx q[2];
rz(-1.3996404) q[2];
sx q[2];
rz(-1.0199176) q[2];
rz(2.6323281) q[3];
sx q[3];
rz(-1.441322) q[3];
sx q[3];
rz(3.063108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7996063) q[0];
sx q[0];
rz(-1.5155563) q[0];
sx q[0];
rz(-1.9827783) q[0];
rz(-0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(1.4656969) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138595) q[0];
sx q[0];
rz(-0.056170551) q[0];
sx q[0];
rz(2.7059731) q[0];
rz(-1.2242555) q[2];
sx q[2];
rz(-0.90991352) q[2];
sx q[2];
rz(-2.5310972) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1701057) q[1];
sx q[1];
rz(-2.4772812) q[1];
sx q[1];
rz(0.12477915) q[1];
x q[2];
rz(1.8578148) q[3];
sx q[3];
rz(-0.94285184) q[3];
sx q[3];
rz(2.7046741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9965685) q[2];
sx q[2];
rz(-1.327876) q[2];
sx q[2];
rz(-0.81929755) q[2];
rz(0.79646349) q[3];
sx q[3];
rz(-0.4070681) q[3];
sx q[3];
rz(1.6527294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(1.2040937) q[0];
sx q[0];
rz(-0.096385328) q[0];
sx q[0];
rz(-0.82164422) q[0];
rz(0.67283982) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(-1.0427262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561036) q[0];
sx q[0];
rz(-0.21767958) q[0];
sx q[0];
rz(-0.34164048) q[0];
rz(-0.94093948) q[2];
sx q[2];
rz(-0.90412882) q[2];
sx q[2];
rz(-2.710944) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4303226) q[1];
sx q[1];
rz(-1.2001598) q[1];
sx q[1];
rz(-2.1881585) q[1];
rz(-pi) q[2];
rz(-0.53570536) q[3];
sx q[3];
rz(-2.8495478) q[3];
sx q[3];
rz(0.10242187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.78160703) q[2];
sx q[2];
rz(-1.2722509) q[2];
sx q[2];
rz(-2.2992415) q[2];
rz(0.53884566) q[3];
sx q[3];
rz(-1.4875965) q[3];
sx q[3];
rz(0.52411383) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6738324) q[0];
sx q[0];
rz(-0.14586511) q[0];
sx q[0];
rz(1.9061331) q[0];
rz(2.1927059) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(-1.6611151) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47280374) q[0];
sx q[0];
rz(-1.2929217) q[0];
sx q[0];
rz(0.86937277) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83809488) q[2];
sx q[2];
rz(-2.3857255) q[2];
sx q[2];
rz(-0.91516337) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64609194) q[1];
sx q[1];
rz(-1.4061855) q[1];
sx q[1];
rz(1.6465523) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32835378) q[3];
sx q[3];
rz(-1.8002276) q[3];
sx q[3];
rz(0.52426527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71198717) q[2];
sx q[2];
rz(-1.4236071) q[2];
sx q[2];
rz(-0.31022662) q[2];
rz(0.45983908) q[3];
sx q[3];
rz(-2.583677) q[3];
sx q[3];
rz(1.7673813) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9643758) q[0];
sx q[0];
rz(-2.0073267) q[0];
sx q[0];
rz(1.0649756) q[0];
rz(-0.89827697) q[1];
sx q[1];
rz(-0.54294642) q[1];
sx q[1];
rz(-0.44566659) q[1];
rz(-0.90625554) q[2];
sx q[2];
rz(-1.8977842) q[2];
sx q[2];
rz(2.1791853) q[2];
rz(-0.38705714) q[3];
sx q[3];
rz(-0.94578065) q[3];
sx q[3];
rz(-2.9972535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
