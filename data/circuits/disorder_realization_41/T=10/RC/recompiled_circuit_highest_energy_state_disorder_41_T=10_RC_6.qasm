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
rz(0.32349411) q[0];
sx q[0];
rz(-0.20322023) q[0];
sx q[0];
rz(0.30015552) q[0];
rz(2.7872941) q[1];
sx q[1];
rz(-1.2343255) q[1];
sx q[1];
rz(-0.77114463) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78438877) q[0];
sx q[0];
rz(-1.3805806) q[0];
sx q[0];
rz(-1.7611544) q[0];
rz(-pi) q[1];
rz(2.6080898) q[2];
sx q[2];
rz(-1.9868104) q[2];
sx q[2];
rz(-0.5678643) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31509545) q[1];
sx q[1];
rz(-2.5719448) q[1];
sx q[1];
rz(0.82514067) q[1];
rz(-pi) q[2];
rz(-3.0740115) q[3];
sx q[3];
rz(-2.0039119) q[3];
sx q[3];
rz(-1.6855296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99042088) q[2];
sx q[2];
rz(-2.6800818) q[2];
sx q[2];
rz(-2.0578461) q[2];
rz(-2.566973) q[3];
sx q[3];
rz(-1.5465522) q[3];
sx q[3];
rz(0.77742022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8271178) q[0];
sx q[0];
rz(-2.5125393) q[0];
sx q[0];
rz(0.3984867) q[0];
rz(-1.1442979) q[1];
sx q[1];
rz(-2.5356348) q[1];
sx q[1];
rz(-0.95432895) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935611) q[0];
sx q[0];
rz(-1.5662417) q[0];
sx q[0];
rz(1.5702374) q[0];
rz(0.91314885) q[2];
sx q[2];
rz(-1.4962915) q[2];
sx q[2];
rz(2.7721921) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0132389) q[1];
sx q[1];
rz(-1.8132105) q[1];
sx q[1];
rz(2.1552312) q[1];
rz(-pi) q[2];
rz(-0.70084533) q[3];
sx q[3];
rz(-0.66674816) q[3];
sx q[3];
rz(-1.4481883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2362471) q[2];
sx q[2];
rz(-0.50062847) q[2];
sx q[2];
rz(-0.11032571) q[2];
rz(-0.8253544) q[3];
sx q[3];
rz(-0.76485157) q[3];
sx q[3];
rz(-2.4013405) q[3];
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
rz(-1.9731307) q[0];
sx q[0];
rz(-2.923712) q[0];
sx q[0];
rz(-2.2595898) q[0];
rz(1.0285671) q[1];
sx q[1];
rz(-2.5879637) q[1];
sx q[1];
rz(0.91241589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9498429) q[0];
sx q[0];
rz(-1.1198567) q[0];
sx q[0];
rz(0.080206932) q[0];
rz(-pi) q[1];
rz(-0.13752748) q[2];
sx q[2];
rz(-2.418387) q[2];
sx q[2];
rz(-2.9971426) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0076651) q[1];
sx q[1];
rz(-1.0093736) q[1];
sx q[1];
rz(2.7473161) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0718143) q[3];
sx q[3];
rz(-2.6423965) q[3];
sx q[3];
rz(0.62455356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2622111) q[2];
sx q[2];
rz(-2.8935581) q[2];
sx q[2];
rz(-0.81825078) q[2];
rz(-2.7815871) q[3];
sx q[3];
rz(-1.8428948) q[3];
sx q[3];
rz(0.037809614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70571947) q[0];
sx q[0];
rz(-2.191045) q[0];
sx q[0];
rz(-1.1608423) q[0];
rz(2.1707161) q[1];
sx q[1];
rz(-2.2258046) q[1];
sx q[1];
rz(-0.11347778) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023983) q[0];
sx q[0];
rz(-1.9736088) q[0];
sx q[0];
rz(2.4871123) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1283896) q[2];
sx q[2];
rz(-1.6062066) q[2];
sx q[2];
rz(-1.7865045) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.27423985) q[1];
sx q[1];
rz(-0.93207658) q[1];
sx q[1];
rz(-2.252905) q[1];
rz(0.7917801) q[3];
sx q[3];
rz(-2.3327475) q[3];
sx q[3];
rz(2.5616982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11161441) q[2];
sx q[2];
rz(-1.3361715) q[2];
sx q[2];
rz(-0.90799904) q[2];
rz(-2.8152483) q[3];
sx q[3];
rz(-0.55485266) q[3];
sx q[3];
rz(-0.42181304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57080764) q[0];
sx q[0];
rz(-0.76199216) q[0];
sx q[0];
rz(1.0170035) q[0];
rz(-0.48655888) q[1];
sx q[1];
rz(-1.0252955) q[1];
sx q[1];
rz(-3.0466381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96387917) q[0];
sx q[0];
rz(-1.6780319) q[0];
sx q[0];
rz(3.043006) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3170065) q[2];
sx q[2];
rz(-1.4207625) q[2];
sx q[2];
rz(2.4730446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50852114) q[1];
sx q[1];
rz(-1.2969368) q[1];
sx q[1];
rz(-0.36063309) q[1];
rz(-pi) q[2];
rz(-1.8848991) q[3];
sx q[3];
rz(-1.8669379) q[3];
sx q[3];
rz(1.2552798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2056556) q[2];
sx q[2];
rz(-0.63465261) q[2];
sx q[2];
rz(0.46979365) q[2];
rz(-0.69822407) q[3];
sx q[3];
rz(-1.2263115) q[3];
sx q[3];
rz(-0.32018143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6093589) q[0];
sx q[0];
rz(-2.4025752) q[0];
sx q[0];
rz(-2.9267689) q[0];
rz(0.23669446) q[1];
sx q[1];
rz(-0.57699811) q[1];
sx q[1];
rz(-2.7223041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474684) q[0];
sx q[0];
rz(-1.6162833) q[0];
sx q[0];
rz(1.4922597) q[0];
x q[1];
rz(2.8307011) q[2];
sx q[2];
rz(-2.0902373) q[2];
sx q[2];
rz(-2.3497557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57945573) q[1];
sx q[1];
rz(-2.750235) q[1];
sx q[1];
rz(-1.0117024) q[1];
x q[2];
rz(0.71281302) q[3];
sx q[3];
rz(-1.4199004) q[3];
sx q[3];
rz(-0.11980443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.145891) q[2];
sx q[2];
rz(-0.13549165) q[2];
sx q[2];
rz(0.091391407) q[2];
rz(-0.35178301) q[3];
sx q[3];
rz(-2.3942949) q[3];
sx q[3];
rz(-2.6763797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2202989) q[0];
sx q[0];
rz(-1.969368) q[0];
sx q[0];
rz(-1.175513) q[0];
rz(0.57918864) q[1];
sx q[1];
rz(-2.3349031) q[1];
sx q[1];
rz(-2.9520292) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2611147) q[0];
sx q[0];
rz(-1.7715329) q[0];
sx q[0];
rz(1.6693382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6734419) q[2];
sx q[2];
rz(-0.70404875) q[2];
sx q[2];
rz(1.2873905) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4259262) q[1];
sx q[1];
rz(-2.3235011) q[1];
sx q[1];
rz(-0.010424213) q[1];
rz(-pi) q[2];
x q[2];
rz(0.060154288) q[3];
sx q[3];
rz(-2.0904448) q[3];
sx q[3];
rz(-1.026767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8976979) q[2];
sx q[2];
rz(-1.2995259) q[2];
sx q[2];
rz(0.50590903) q[2];
rz(-0.77887744) q[3];
sx q[3];
rz(-2.7454822) q[3];
sx q[3];
rz(1.8469384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8706354) q[0];
sx q[0];
rz(-2.0896235) q[0];
sx q[0];
rz(-1.9860995) q[0];
rz(-0.0011860154) q[1];
sx q[1];
rz(-0.78758925) q[1];
sx q[1];
rz(-0.40518618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4219932) q[0];
sx q[0];
rz(-2.036839) q[0];
sx q[0];
rz(-2.5525722) q[0];
x q[1];
rz(1.979902) q[2];
sx q[2];
rz(-1.7120581) q[2];
sx q[2];
rz(-2.2337332) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2470513) q[1];
sx q[1];
rz(-0.76746583) q[1];
sx q[1];
rz(1.6805499) q[1];
x q[2];
rz(-3.013582) q[3];
sx q[3];
rz(-1.6388005) q[3];
sx q[3];
rz(-1.3101414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.323641) q[2];
sx q[2];
rz(-1.9743071) q[2];
sx q[2];
rz(-0.52158588) q[2];
rz(-3.0449872) q[3];
sx q[3];
rz(-2.1229027) q[3];
sx q[3];
rz(2.6507613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38799649) q[0];
sx q[0];
rz(-2.2612408) q[0];
sx q[0];
rz(-3.0269347) q[0];
rz(-1.313734) q[1];
sx q[1];
rz(-1.6856245) q[1];
sx q[1];
rz(-0.1350666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60914441) q[0];
sx q[0];
rz(-1.3125064) q[0];
sx q[0];
rz(1.7946971) q[0];
rz(-1.3038605) q[2];
sx q[2];
rz(-1.1914754) q[2];
sx q[2];
rz(-2.1960559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2429906) q[1];
sx q[1];
rz(-0.95546104) q[1];
sx q[1];
rz(-0.19886417) q[1];
rz(2.1296639) q[3];
sx q[3];
rz(-1.5084576) q[3];
sx q[3];
rz(-1.896794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8257398) q[2];
sx q[2];
rz(-1.9261381) q[2];
sx q[2];
rz(-2.1960171) q[2];
rz(0.25997508) q[3];
sx q[3];
rz(-2.1187449) q[3];
sx q[3];
rz(-2.0135349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0876227) q[0];
sx q[0];
rz(-2.051351) q[0];
sx q[0];
rz(-1.7386275) q[0];
rz(-2.2258017) q[1];
sx q[1];
rz(-0.61640888) q[1];
sx q[1];
rz(2.2260407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957324) q[0];
sx q[0];
rz(-1.440319) q[0];
sx q[0];
rz(0.097921485) q[0];
rz(-1.8981839) q[2];
sx q[2];
rz(-0.54556525) q[2];
sx q[2];
rz(2.4107673) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7448235) q[1];
sx q[1];
rz(-1.9717378) q[1];
sx q[1];
rz(2.045782) q[1];
rz(-2.2804349) q[3];
sx q[3];
rz(-1.7000515) q[3];
sx q[3];
rz(-1.2354163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.12330684) q[2];
sx q[2];
rz(-2.4985963) q[2];
sx q[2];
rz(0.40261889) q[2];
rz(-1.9649402) q[3];
sx q[3];
rz(-0.28665701) q[3];
sx q[3];
rz(0.76796842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1422414) q[0];
sx q[0];
rz(-2.169123) q[0];
sx q[0];
rz(2.1219667) q[0];
rz(0.6877407) q[1];
sx q[1];
rz(-0.80324695) q[1];
sx q[1];
rz(-1.3245503) q[1];
rz(2.3024846) q[2];
sx q[2];
rz(-1.0232836) q[2];
sx q[2];
rz(-2.8689775) q[2];
rz(0.34360828) q[3];
sx q[3];
rz(-0.6597198) q[3];
sx q[3];
rz(2.3312286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
