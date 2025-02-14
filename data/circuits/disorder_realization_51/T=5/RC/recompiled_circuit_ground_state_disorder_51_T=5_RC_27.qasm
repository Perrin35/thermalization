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
rz(-0.0039761877) q[0];
sx q[0];
rz(-3.0206326) q[0];
rz(0.51789415) q[1];
sx q[1];
rz(2.8105812) q[1];
sx q[1];
rz(10.612648) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5802104) q[0];
sx q[0];
rz(-0.31162173) q[0];
sx q[0];
rz(-2.8557957) q[0];
rz(-pi) q[1];
rz(2.0969432) q[2];
sx q[2];
rz(-1.4441212) q[2];
sx q[2];
rz(1.0110477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.647794) q[1];
sx q[1];
rz(-1.0435074) q[1];
sx q[1];
rz(2.6721558) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.066066001) q[3];
sx q[3];
rz(-0.63263901) q[3];
sx q[3];
rz(0.047601117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0821705) q[2];
sx q[2];
rz(-0.94258451) q[2];
sx q[2];
rz(-1.7131294) q[2];
rz(1.7510471) q[3];
sx q[3];
rz(-2.3640552) q[3];
sx q[3];
rz(-0.20974222) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9222337) q[0];
sx q[0];
rz(-1.5809504) q[0];
sx q[0];
rz(-1.6844164) q[0];
rz(0.68151418) q[1];
sx q[1];
rz(-0.69029713) q[1];
sx q[1];
rz(-0.95726475) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831612) q[0];
sx q[0];
rz(-0.51379075) q[0];
sx q[0];
rz(-2.7130037) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.817746) q[2];
sx q[2];
rz(-1.3829573) q[2];
sx q[2];
rz(-2.2855594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71979501) q[1];
sx q[1];
rz(-1.6198427) q[1];
sx q[1];
rz(-1.1625467) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010250895) q[3];
sx q[3];
rz(-1.5019377) q[3];
sx q[3];
rz(-2.0673925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7822781) q[2];
sx q[2];
rz(-2.094153) q[2];
sx q[2];
rz(0.90685833) q[2];
rz(-2.4685229) q[3];
sx q[3];
rz(-2.7570717) q[3];
sx q[3];
rz(-1.8460021) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68536711) q[0];
sx q[0];
rz(-0.65995589) q[0];
sx q[0];
rz(0.54687706) q[0];
rz(1.7710255) q[1];
sx q[1];
rz(-1.9799045) q[1];
sx q[1];
rz(0.10017698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91104394) q[0];
sx q[0];
rz(-1.3063138) q[0];
sx q[0];
rz(3.0039211) q[0];
x q[1];
rz(-0.64073784) q[2];
sx q[2];
rz(-0.60534873) q[2];
sx q[2];
rz(-1.5126765) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2168536) q[1];
sx q[1];
rz(-1.3432055) q[1];
sx q[1];
rz(0.74149577) q[1];
rz(1.1421575) q[3];
sx q[3];
rz(-1.1446307) q[3];
sx q[3];
rz(-1.4166735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3935516) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(2.098341) q[2];
rz(-0.68759632) q[3];
sx q[3];
rz(-0.50191534) q[3];
sx q[3];
rz(2.8858394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2954751) q[0];
sx q[0];
rz(-0.44545528) q[0];
sx q[0];
rz(-2.1050146) q[0];
rz(-1.8702033) q[1];
sx q[1];
rz(-0.82415736) q[1];
sx q[1];
rz(-1.8545256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6223479) q[0];
sx q[0];
rz(-2.17832) q[0];
sx q[0];
rz(0.95444478) q[0];
rz(2.4802568) q[2];
sx q[2];
rz(-1.9171538) q[2];
sx q[2];
rz(-0.48717425) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92136473) q[1];
sx q[1];
rz(-1.2829687) q[1];
sx q[1];
rz(-2.1975321) q[1];
rz(-pi) q[2];
rz(0.98388328) q[3];
sx q[3];
rz(-0.5183087) q[3];
sx q[3];
rz(-0.89215088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6048772) q[2];
sx q[2];
rz(-0.36467364) q[2];
sx q[2];
rz(-3.1216915) q[2];
rz(2.5455635) q[3];
sx q[3];
rz(-2.8607131) q[3];
sx q[3];
rz(-1.9635487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48536456) q[0];
sx q[0];
rz(-1.8940268) q[0];
sx q[0];
rz(1.1744936) q[0];
rz(-0.30612048) q[1];
sx q[1];
rz(-1.0447634) q[1];
sx q[1];
rz(-1.7721446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58254439) q[0];
sx q[0];
rz(-1.567548) q[0];
sx q[0];
rz(1.5772468) q[0];
rz(-pi) q[1];
rz(-2.1691983) q[2];
sx q[2];
rz(-1.7142503) q[2];
sx q[2];
rz(-0.17696887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9766531) q[1];
sx q[1];
rz(-2.7639228) q[1];
sx q[1];
rz(-1.151848) q[1];
rz(-pi) q[2];
rz(-1.1569275) q[3];
sx q[3];
rz(-2.0479432) q[3];
sx q[3];
rz(-2.8780495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8451763) q[2];
sx q[2];
rz(-2.3974819) q[2];
sx q[2];
rz(2.3386492) q[2];
rz(1.1154277) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(1.7723134) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.707412) q[0];
sx q[0];
rz(-0.65874225) q[0];
sx q[0];
rz(-2.1408268) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66866691) q[2];
sx q[2];
rz(-1.200182) q[2];
sx q[2];
rz(0.47793717) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2032092) q[1];
sx q[1];
rz(-0.46493772) q[1];
sx q[1];
rz(0.26588456) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93490113) q[3];
sx q[3];
rz(-0.29014698) q[3];
sx q[3];
rz(-1.6573418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4623146) q[2];
sx q[2];
rz(-2.9913112) q[2];
sx q[2];
rz(1.3432937) q[2];
rz(0.51370931) q[3];
sx q[3];
rz(-1.6535583) q[3];
sx q[3];
rz(0.59521365) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3419471) q[0];
sx q[0];
rz(-1.8193614) q[0];
sx q[0];
rz(-0.42022002) q[0];
rz(1.7814024) q[1];
sx q[1];
rz(-1.1667292) q[1];
sx q[1];
rz(0.80418783) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5468381) q[0];
sx q[0];
rz(-1.9512984) q[0];
sx q[0];
rz(0.98989931) q[0];
rz(-pi) q[1];
rz(-1.8558414) q[2];
sx q[2];
rz(-2.3346296) q[2];
sx q[2];
rz(-2.3318044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7090164) q[1];
sx q[1];
rz(-2.1423999) q[1];
sx q[1];
rz(-2.9663621) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7358094) q[3];
sx q[3];
rz(-1.6783774) q[3];
sx q[3];
rz(3.0522814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5899902) q[2];
sx q[2];
rz(-1.0162105) q[2];
sx q[2];
rz(1.1691947) q[2];
rz(-0.17360887) q[3];
sx q[3];
rz(-2.2326525) q[3];
sx q[3];
rz(-1.9420697) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46840295) q[0];
sx q[0];
rz(-1.0262187) q[0];
sx q[0];
rz(1.9386559) q[0];
rz(-2.9007593) q[1];
sx q[1];
rz(-0.92643654) q[1];
sx q[1];
rz(-2.208272) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881977) q[0];
sx q[0];
rz(-1.9472593) q[0];
sx q[0];
rz(0.36247902) q[0];
rz(-pi) q[1];
rz(-2.1432518) q[2];
sx q[2];
rz(-1.7939374) q[2];
sx q[2];
rz(2.9024692) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5211903) q[1];
sx q[1];
rz(-0.4333843) q[1];
sx q[1];
rz(-2.6492725) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15214129) q[3];
sx q[3];
rz(-2.1295432) q[3];
sx q[3];
rz(1.5882176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0241218) q[2];
sx q[2];
rz(-1.6880219) q[2];
sx q[2];
rz(-3.0061099) q[2];
rz(2.1425715) q[3];
sx q[3];
rz(-2.8935367) q[3];
sx q[3];
rz(-2.5854056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67021543) q[0];
sx q[0];
rz(-1.1286292) q[0];
sx q[0];
rz(-1.7673329) q[0];
rz(-2.0595835) q[1];
sx q[1];
rz(-2.3550985) q[1];
sx q[1];
rz(1.5622004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53331965) q[0];
sx q[0];
rz(-1.0315622) q[0];
sx q[0];
rz(3.1031552) q[0];
x q[1];
rz(-0.5600881) q[2];
sx q[2];
rz(-1.9011407) q[2];
sx q[2];
rz(2.7719088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1154707) q[1];
sx q[1];
rz(-1.9023832) q[1];
sx q[1];
rz(-2.0381171) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0252941) q[3];
sx q[3];
rz(-0.70279944) q[3];
sx q[3];
rz(-1.1490319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18206) q[2];
sx q[2];
rz(-1.5740732) q[2];
sx q[2];
rz(0.62425295) q[2];
rz(-2.4652081) q[3];
sx q[3];
rz(-1.1712149) q[3];
sx q[3];
rz(-0.83522767) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71037173) q[0];
sx q[0];
rz(-1.2483163) q[0];
sx q[0];
rz(-1.3488019) q[0];
rz(2.8377332) q[1];
sx q[1];
rz(-1.6108578) q[1];
sx q[1];
rz(0.87711984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9187885) q[0];
sx q[0];
rz(-0.74367803) q[0];
sx q[0];
rz(-0.2144465) q[0];
rz(1.2459846) q[2];
sx q[2];
rz(-1.6742252) q[2];
sx q[2];
rz(1.499044) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34906975) q[1];
sx q[1];
rz(-1.190509) q[1];
sx q[1];
rz(-1.1009586) q[1];
x q[2];
rz(0.95939674) q[3];
sx q[3];
rz(-1.093321) q[3];
sx q[3];
rz(1.9408664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0861686) q[2];
sx q[2];
rz(-0.40579) q[2];
sx q[2];
rz(1.0322734) q[2];
rz(-2.0530733) q[3];
sx q[3];
rz(-1.4884596) q[3];
sx q[3];
rz(-2.3672339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9536204) q[0];
sx q[0];
rz(-0.80360501) q[0];
sx q[0];
rz(-3.0927717) q[0];
rz(-1.8232518) q[1];
sx q[1];
rz(-1.4069469) q[1];
sx q[1];
rz(-2.5896752) q[1];
rz(0.55808347) q[2];
sx q[2];
rz(-0.48102745) q[2];
sx q[2];
rz(-0.68651909) q[2];
rz(-1.134675) q[3];
sx q[3];
rz(-1.5584617) q[3];
sx q[3];
rz(-0.8212318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
