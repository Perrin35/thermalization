OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(-5.6145515) q[1];
sx q[1];
rz(0.86548391) q[1];
sx q[1];
rz(15.794985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4092769) q[0];
sx q[0];
rz(-1.1999994) q[0];
sx q[0];
rz(-0.91863527) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1331698) q[2];
sx q[2];
rz(-0.77169092) q[2];
sx q[2];
rz(2.0193677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0793003) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(-2.9855707) q[1];
rz(-pi) q[2];
rz(3.0356785) q[3];
sx q[3];
rz(-0.4142136) q[3];
sx q[3];
rz(-0.29990444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0027851) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.194838) q[0];
rz(-3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502991) q[0];
sx q[0];
rz(-2.5851558) q[0];
sx q[0];
rz(0.87470212) q[0];
x q[1];
rz(1.8325009) q[2];
sx q[2];
rz(-1.3771025) q[2];
sx q[2];
rz(1.9270093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.61987585) q[1];
sx q[1];
rz(-2.4881425) q[1];
sx q[1];
rz(2.7481649) q[1];
rz(-0.46044465) q[3];
sx q[3];
rz(-1.5213955) q[3];
sx q[3];
rz(-1.3747665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(2.3475032) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2394543) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(-2.7368271) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.3084897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.20298) q[0];
sx q[0];
rz(-0.70883195) q[0];
sx q[0];
rz(0.75074408) q[0];
rz(-pi) q[1];
rz(-0.38396118) q[2];
sx q[2];
rz(-2.5703891) q[2];
sx q[2];
rz(0.9108033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82061003) q[1];
sx q[1];
rz(-2.2690986) q[1];
sx q[1];
rz(0.33336158) q[1];
rz(-2.7029196) q[3];
sx q[3];
rz(-1.6524501) q[3];
sx q[3];
rz(-2.4077533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.219316) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(-2.0698047) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(-2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(2.2256057) q[0];
rz(-2.6782716) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076544882) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(1.5426427) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3154503) q[2];
sx q[2];
rz(-2.38378) q[2];
sx q[2];
rz(0.10105029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.108236) q[1];
sx q[1];
rz(-1.557096) q[1];
sx q[1];
rz(-1.800888) q[1];
x q[2];
rz(-1.7980868) q[3];
sx q[3];
rz(-2.3152581) q[3];
sx q[3];
rz(0.28971653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89381924) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(-1.4438859) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(-1.7516288) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(-0.46868971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5422573) q[0];
sx q[0];
rz(-2.242803) q[0];
sx q[0];
rz(0.50962781) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52207559) q[2];
sx q[2];
rz(-2.7895088) q[2];
sx q[2];
rz(-3.0662231) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3846181) q[1];
sx q[1];
rz(-2.9006835) q[1];
sx q[1];
rz(0.5730281) q[1];
rz(-0.080805578) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9995352) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-0.86432499) q[2];
rz(-0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(-2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(0.19700225) q[0];
rz(-1.7794094) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(0.33624712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42423781) q[0];
sx q[0];
rz(-0.80168085) q[0];
sx q[0];
rz(2.5108811) q[0];
rz(-pi) q[1];
rz(0.59351144) q[2];
sx q[2];
rz(-2.0715908) q[2];
sx q[2];
rz(-1.8261432) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5583401) q[1];
sx q[1];
rz(-1.5234158) q[1];
sx q[1];
rz(-2.5517795) q[1];
x q[2];
rz(-0.20506649) q[3];
sx q[3];
rz(-0.70687095) q[3];
sx q[3];
rz(-0.52458483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53529915) q[2];
sx q[2];
rz(-1.665411) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(-1.2396631) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(-0.25892648) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(2.0475725) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0897113) q[0];
sx q[0];
rz(-1.9117172) q[0];
sx q[0];
rz(-1.107035) q[0];
rz(-pi) q[1];
x q[1];
rz(0.012374087) q[2];
sx q[2];
rz(-1.7933328) q[2];
sx q[2];
rz(-2.6518133) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3257521) q[1];
sx q[1];
rz(-0.098512352) q[1];
sx q[1];
rz(1.2835531) q[1];
rz(-pi) q[2];
rz(-0.098071531) q[3];
sx q[3];
rz(-1.5925928) q[3];
sx q[3];
rz(1.9779713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(-0.32067498) q[2];
rz(2.705412) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(2.5514305) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(2.337713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76669508) q[0];
sx q[0];
rz(-0.98537966) q[0];
sx q[0];
rz(1.5234408) q[0];
rz(1.3857533) q[2];
sx q[2];
rz(-2.1171326) q[2];
sx q[2];
rz(1.1748558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4623973) q[1];
sx q[1];
rz(-0.089226626) q[1];
sx q[1];
rz(-1.6855082) q[1];
rz(-pi) q[2];
rz(1.4804833) q[3];
sx q[3];
rz(-2.2323425) q[3];
sx q[3];
rz(-0.46935287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-1.0104898) q[2];
rz(-0.60339749) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-0.4075152) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(0.75072748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84839612) q[0];
sx q[0];
rz(-1.2788532) q[0];
sx q[0];
rz(-1.8672724) q[0];
rz(-pi) q[1];
rz(-1.4892011) q[2];
sx q[2];
rz(-1.2921234) q[2];
sx q[2];
rz(2.1747053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6299906) q[1];
sx q[1];
rz(-1.2641608) q[1];
sx q[1];
rz(0.2125477) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.084089355) q[3];
sx q[3];
rz(-0.30232271) q[3];
sx q[3];
rz(1.38894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0687381) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(-1.1661952) q[2];
rz(-1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5831379) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(-1.7512084) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(1.6814544) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64797348) q[0];
sx q[0];
rz(-2.4612392) q[0];
sx q[0];
rz(2.2749659) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5286469) q[2];
sx q[2];
rz(-0.50439207) q[2];
sx q[2];
rz(1.1500037) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8148988) q[1];
sx q[1];
rz(-1.9683451) q[1];
sx q[1];
rz(-1.7705998) q[1];
rz(-pi) q[2];
x q[2];
rz(2.371993) q[3];
sx q[3];
rz(-1.0255314) q[3];
sx q[3];
rz(-1.256497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(-1.6798518) q[2];
rz(-1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5158952) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(-1.3810146) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(-1.9258326) q[2];
sx q[2];
rz(-2.3264865) q[2];
sx q[2];
rz(-2.9441499) q[2];
rz(-0.67301987) q[3];
sx q[3];
rz(-1.2613847) q[3];
sx q[3];
rz(-3.1223084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];