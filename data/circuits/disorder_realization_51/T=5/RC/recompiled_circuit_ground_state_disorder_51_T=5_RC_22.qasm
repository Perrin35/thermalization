OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.12282523) q[0];
sx q[0];
rz(3.1376165) q[0];
sx q[0];
rz(9.3038179) q[0];
rz(0.51789415) q[1];
sx q[1];
rz(-0.33101141) q[1];
sx q[1];
rz(1.9537227) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56138229) q[0];
sx q[0];
rz(-0.31162173) q[0];
sx q[0];
rz(0.28579692) q[0];
x q[1];
rz(1.0446494) q[2];
sx q[2];
rz(-1.4441212) q[2];
sx q[2];
rz(2.130545) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49379865) q[1];
sx q[1];
rz(-2.0980853) q[1];
sx q[1];
rz(2.6721558) q[1];
x q[2];
rz(0.066066001) q[3];
sx q[3];
rz(-0.63263901) q[3];
sx q[3];
rz(3.0939915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0594222) q[2];
sx q[2];
rz(-2.1990081) q[2];
sx q[2];
rz(-1.4284632) q[2];
rz(-1.7510471) q[3];
sx q[3];
rz(-0.7775375) q[3];
sx q[3];
rz(2.9318504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9222337) q[0];
sx q[0];
rz(-1.5809504) q[0];
sx q[0];
rz(1.4571762) q[0];
rz(-0.68151418) q[1];
sx q[1];
rz(-2.4512955) q[1];
sx q[1];
rz(-0.95726475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5663366) q[0];
sx q[0];
rz(-1.1073551) q[0];
sx q[0];
rz(1.3404247) q[0];
rz(-pi) q[1];
x q[1];
rz(0.817746) q[2];
sx q[2];
rz(-1.7586353) q[2];
sx q[2];
rz(0.85603324) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4217976) q[1];
sx q[1];
rz(-1.52175) q[1];
sx q[1];
rz(-1.979046) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1313418) q[3];
sx q[3];
rz(-1.5019377) q[3];
sx q[3];
rz(1.0742002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35931453) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(-0.90685833) q[2];
rz(-2.4685229) q[3];
sx q[3];
rz(-0.38452092) q[3];
sx q[3];
rz(1.8460021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68536711) q[0];
sx q[0];
rz(-2.4816368) q[0];
sx q[0];
rz(-2.5947156) q[0];
rz(-1.3705672) q[1];
sx q[1];
rz(-1.1616881) q[1];
sx q[1];
rz(-0.10017698) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2305487) q[0];
sx q[0];
rz(-1.3063138) q[0];
sx q[0];
rz(-3.0039211) q[0];
x q[1];
rz(-2.6351026) q[2];
sx q[2];
rz(-1.9178941) q[2];
sx q[2];
rz(-2.533415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.924739) q[1];
sx q[1];
rz(-1.3432055) q[1];
sx q[1];
rz(-0.74149577) q[1];
rz(-1.1421575) q[3];
sx q[3];
rz(-1.1446307) q[3];
sx q[3];
rz(-1.7249191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.748041) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(-2.098341) q[2];
rz(-2.4539963) q[3];
sx q[3];
rz(-0.50191534) q[3];
sx q[3];
rz(0.25575328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8461175) q[0];
sx q[0];
rz(-2.6961374) q[0];
sx q[0];
rz(-1.0365781) q[0];
rz(1.2713894) q[1];
sx q[1];
rz(-2.3174353) q[1];
sx q[1];
rz(1.8545256) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7298147) q[0];
sx q[0];
rz(-2.3049816) q[0];
sx q[0];
rz(-2.4479624) q[0];
rz(-pi) q[1];
rz(-0.53127473) q[2];
sx q[2];
rz(-2.4073045) q[2];
sx q[2];
rz(-1.6466494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8657664) q[1];
sx q[1];
rz(-2.4600907) q[1];
sx q[1];
rz(1.1033415) q[1];
rz(2.014145) q[3];
sx q[3];
rz(-1.2928767) q[3];
sx q[3];
rz(-1.2025361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6048772) q[2];
sx q[2];
rz(-2.776919) q[2];
sx q[2];
rz(0.019901179) q[2];
rz(0.59602916) q[3];
sx q[3];
rz(-0.28087956) q[3];
sx q[3];
rz(1.178044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48536456) q[0];
sx q[0];
rz(-1.2475659) q[0];
sx q[0];
rz(-1.9670991) q[0];
rz(2.8354722) q[1];
sx q[1];
rz(-1.0447634) q[1];
sx q[1];
rz(1.3694481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1533617) q[0];
sx q[0];
rz(-1.5643459) q[0];
sx q[0];
rz(3.1383443) q[0];
rz(0.97239437) q[2];
sx q[2];
rz(-1.4273424) q[2];
sx q[2];
rz(-2.9646238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3432797) q[1];
sx q[1];
rz(-1.4202182) q[1];
sx q[1];
rz(1.9184789) q[1];
rz(-1.1569275) q[3];
sx q[3];
rz(-2.0479432) q[3];
sx q[3];
rz(0.26354313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8451763) q[2];
sx q[2];
rz(-0.74411074) q[2];
sx q[2];
rz(-0.80294341) q[2];
rz(-2.0261649) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(1.7723134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1882316) q[0];
sx q[0];
rz(-0.53934923) q[0];
sx q[0];
rz(-0.11235919) q[0];
rz(-0.27680963) q[1];
sx q[1];
rz(-1.9748297) q[1];
sx q[1];
rz(-2.4940431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.707412) q[0];
sx q[0];
rz(-0.65874225) q[0];
sx q[0];
rz(1.0007658) q[0];
rz(2.0306088) q[2];
sx q[2];
rz(-0.95470482) q[2];
sx q[2];
rz(1.3715708) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6426443) q[1];
sx q[1];
rz(-2.0181839) q[1];
sx q[1];
rz(1.7018464) q[1];
rz(-pi) q[2];
rz(-2.2066915) q[3];
sx q[3];
rz(-0.29014698) q[3];
sx q[3];
rz(1.4842509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4623146) q[2];
sx q[2];
rz(-2.9913112) q[2];
sx q[2];
rz(-1.798299) q[2];
rz(-2.6278833) q[3];
sx q[3];
rz(-1.6535583) q[3];
sx q[3];
rz(-2.546379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3419471) q[0];
sx q[0];
rz(-1.3222313) q[0];
sx q[0];
rz(-2.7213726) q[0];
rz(-1.3601903) q[1];
sx q[1];
rz(-1.9748634) q[1];
sx q[1];
rz(-0.80418783) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73690945) q[0];
sx q[0];
rz(-1.0361639) q[0];
sx q[0];
rz(2.695309) q[0];
rz(-pi) q[1];
rz(-1.8558414) q[2];
sx q[2];
rz(-2.3346296) q[2];
sx q[2];
rz(-2.3318044) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11629352) q[1];
sx q[1];
rz(-2.5466047) q[1];
sx q[1];
rz(1.8354675) q[1];
x q[2];
rz(0.10905091) q[3];
sx q[3];
rz(-1.4067459) q[3];
sx q[3];
rz(-1.6422288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5516025) q[2];
sx q[2];
rz(-2.1253822) q[2];
sx q[2];
rz(-1.1691947) q[2];
rz(-0.17360887) q[3];
sx q[3];
rz(-0.9089402) q[3];
sx q[3];
rz(-1.199523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6731897) q[0];
sx q[0];
rz(-1.0262187) q[0];
sx q[0];
rz(-1.9386559) q[0];
rz(2.9007593) q[1];
sx q[1];
rz(-0.92643654) q[1];
sx q[1];
rz(-0.9333207) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6872662) q[0];
sx q[0];
rz(-2.6250702) q[0];
sx q[0];
rz(0.83965001) q[0];
rz(-0.26367311) q[2];
sx q[2];
rz(-2.1273551) q[2];
sx q[2];
rz(1.9515748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1541862) q[1];
sx q[1];
rz(-1.1917104) q[1];
sx q[1];
rz(-1.7861219) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8086173) q[3];
sx q[3];
rz(-0.57695848) q[3];
sx q[3];
rz(1.2718309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11747083) q[2];
sx q[2];
rz(-1.4535707) q[2];
sx q[2];
rz(3.0061099) q[2];
rz(-2.1425715) q[3];
sx q[3];
rz(-0.24805598) q[3];
sx q[3];
rz(0.55618709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67021543) q[0];
sx q[0];
rz(-1.1286292) q[0];
sx q[0];
rz(1.7673329) q[0];
rz(2.0595835) q[1];
sx q[1];
rz(-0.78649414) q[1];
sx q[1];
rz(1.5622004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0572206) q[0];
sx q[0];
rz(-1.5378152) q[0];
sx q[0];
rz(1.0312366) q[0];
rz(2.5815046) q[2];
sx q[2];
rz(-1.240452) q[2];
sx q[2];
rz(0.3696839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97171451) q[1];
sx q[1];
rz(-2.5757667) q[1];
sx q[1];
rz(-2.2233655) q[1];
x q[2];
rz(1.1162986) q[3];
sx q[3];
rz(-0.70279944) q[3];
sx q[3];
rz(1.9925607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9595327) q[2];
sx q[2];
rz(-1.5740732) q[2];
sx q[2];
rz(-0.62425295) q[2];
rz(0.67638451) q[3];
sx q[3];
rz(-1.1712149) q[3];
sx q[3];
rz(2.306365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312209) q[0];
sx q[0];
rz(-1.2483163) q[0];
sx q[0];
rz(-1.7927908) q[0];
rz(0.30385941) q[1];
sx q[1];
rz(-1.6108578) q[1];
sx q[1];
rz(2.2644728) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064909272) q[0];
sx q[0];
rz(-0.8479894) q[0];
sx q[0];
rz(-1.377489) q[0];
rz(-pi) q[1];
rz(-0.10909144) q[2];
sx q[2];
rz(-1.8938091) q[2];
sx q[2];
rz(-3.035088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5511131) q[1];
sx q[1];
rz(-2.5462228) q[1];
sx q[1];
rz(-2.2941085) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1821959) q[3];
sx q[3];
rz(-1.093321) q[3];
sx q[3];
rz(1.9408664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0861686) q[2];
sx q[2];
rz(-2.7358027) q[2];
sx q[2];
rz(-2.1093192) q[2];
rz(-2.0530733) q[3];
sx q[3];
rz(-1.4884596) q[3];
sx q[3];
rz(0.77435875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9536204) q[0];
sx q[0];
rz(-2.3379876) q[0];
sx q[0];
rz(0.048820989) q[0];
rz(1.3183409) q[1];
sx q[1];
rz(-1.4069469) q[1];
sx q[1];
rz(-2.5896752) q[1];
rz(-1.3011408) q[2];
sx q[2];
rz(-1.9741304) q[2];
sx q[2];
rz(-1.3000549) q[2];
rz(-1.5999888) q[3];
sx q[3];
rz(-2.7053082) q[3];
sx q[3];
rz(-2.3655688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
