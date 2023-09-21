OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(-2.2111501) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(-2.0071109) q[1];
sx q[1];
rz(-1.1073444) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005172) q[0];
sx q[0];
rz(-2.374875) q[0];
sx q[0];
rz(-0.47857743) q[0];
rz(-pi) q[1];
rz(1.0385752) q[2];
sx q[2];
rz(-2.3468809) q[2];
sx q[2];
rz(2.6796535) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2335637) q[1];
sx q[1];
rz(-1.4587914) q[1];
sx q[1];
rz(-1.2909375) q[1];
rz(-pi) q[2];
rz(-2.2060478) q[3];
sx q[3];
rz(-1.8555292) q[3];
sx q[3];
rz(2.2176544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.063623108) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.3280274) q[2];
rz(-2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-0.13197556) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(1.2940787) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(-2.3243288) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0016804455) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(1.2225479) q[0];
x q[1];
rz(-2.0613282) q[2];
sx q[2];
rz(-1.763952) q[2];
sx q[2];
rz(-0.9888538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0664132) q[1];
sx q[1];
rz(-1.315409) q[1];
sx q[1];
rz(-0.65402072) q[1];
rz(1.9278139) q[3];
sx q[3];
rz(-0.5001874) q[3];
sx q[3];
rz(-0.57750765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91784224) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(2.7775653) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(-1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.7746975) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(-2.1752775) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(-1.4470709) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9896302) q[0];
sx q[0];
rz(-1.785779) q[0];
sx q[0];
rz(-1.5887898) q[0];
x q[1];
rz(1.3005199) q[2];
sx q[2];
rz(-2.8405361) q[2];
sx q[2];
rz(-1.9384055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9008873) q[1];
sx q[1];
rz(-1.1966238) q[1];
sx q[1];
rz(0.64970533) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8964642) q[3];
sx q[3];
rz(-1.9760625) q[3];
sx q[3];
rz(0.95642904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7029999) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(-1.1052216) q[2];
rz(-2.3953719) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063909) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(-1.6756469) q[0];
rz(0.28494596) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(2.7526061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3272414) q[0];
sx q[0];
rz(-1.2885433) q[0];
sx q[0];
rz(-1.3319356) q[0];
rz(2.94611) q[2];
sx q[2];
rz(-1.0921548) q[2];
sx q[2];
rz(0.10749707) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4219141) q[1];
sx q[1];
rz(-0.72243566) q[1];
sx q[1];
rz(2.4721485) q[1];
rz(2.9205434) q[3];
sx q[3];
rz(-2.3813558) q[3];
sx q[3];
rz(-1.5628712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.1636461) q[2];
sx q[2];
rz(0.45670613) q[2];
rz(-1.5151954) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(-2.8592498) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.221955) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(-2.4941764) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(2.6470851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9279328) q[0];
sx q[0];
rz(-1.4466009) q[0];
sx q[0];
rz(0.28161063) q[0];
rz(1.7062543) q[2];
sx q[2];
rz(-1.5317305) q[2];
sx q[2];
rz(2.1087697) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33489409) q[1];
sx q[1];
rz(-0.94167275) q[1];
sx q[1];
rz(1.4966399) q[1];
rz(2.6974929) q[3];
sx q[3];
rz(-1.5937623) q[3];
sx q[3];
rz(-0.75957739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(-2.8430856) q[2];
rz(0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(-2.503094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1453778) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(0.19590713) q[0];
rz(0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(-1.235199) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053947833) q[0];
sx q[0];
rz(-1.2637648) q[0];
sx q[0];
rz(-0.35015492) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9094798) q[2];
sx q[2];
rz(-0.89502305) q[2];
sx q[2];
rz(2.0896926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0732121) q[1];
sx q[1];
rz(-2.3975323) q[1];
sx q[1];
rz(0.63853227) q[1];
x q[2];
rz(2.9052832) q[3];
sx q[3];
rz(-1.6522539) q[3];
sx q[3];
rz(2.8636275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6482676) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(-1.7290944) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(0.13599642) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(0.21513367) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(1.1134061) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68590251) q[0];
sx q[0];
rz(-2.6705574) q[0];
sx q[0];
rz(2.5409565) q[0];
x q[1];
rz(1.8298803) q[2];
sx q[2];
rz(-1.2652745) q[2];
sx q[2];
rz(2.6079026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6326633) q[1];
sx q[1];
rz(-0.76970657) q[1];
sx q[1];
rz(2.6973004) q[1];
rz(-1.1235085) q[3];
sx q[3];
rz(-1.7637858) q[3];
sx q[3];
rz(2.6393294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1404184) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(-0.12602885) q[2];
rz(-1.0472939) q[3];
sx q[3];
rz(-1.8208561) q[3];
sx q[3];
rz(-0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814608) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(-1.8966282) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(1.9326899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93223244) q[0];
sx q[0];
rz(-2.8673842) q[0];
sx q[0];
rz(1.8673531) q[0];
rz(-2.2751341) q[2];
sx q[2];
rz(-2.2668215) q[2];
sx q[2];
rz(-2.8528086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.804467) q[1];
sx q[1];
rz(-0.85782385) q[1];
sx q[1];
rz(-1.1285524) q[1];
x q[2];
rz(1.2743837) q[3];
sx q[3];
rz(-2.489438) q[3];
sx q[3];
rz(1.46331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(1.3195999) q[2];
rz(-2.54946) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(-0.035141703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9086583) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(1.3611025) q[0];
rz(1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-0.39168721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4995394) q[0];
sx q[0];
rz(-1.5766489) q[0];
sx q[0];
rz(0.43453479) q[0];
rz(-2.6892745) q[2];
sx q[2];
rz(-1.7512133) q[2];
sx q[2];
rz(-2.7197321) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0605676) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(0.4391567) q[1];
rz(-pi) q[2];
rz(-3.016032) q[3];
sx q[3];
rz(-0.52913044) q[3];
sx q[3];
rz(2.373113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52035511) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(1.3367782) q[2];
rz(-2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(0.80037642) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8005463) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(2.9796519) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25625944) q[0];
sx q[0];
rz(-2.4618039) q[0];
sx q[0];
rz(-1.1023561) q[0];
rz(1.3822046) q[2];
sx q[2];
rz(-1.57975) q[2];
sx q[2];
rz(2.9262275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83878126) q[1];
sx q[1];
rz(-1.4421842) q[1];
sx q[1];
rz(1.7437115) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5374244) q[3];
sx q[3];
rz(-1.0479234) q[3];
sx q[3];
rz(3.0748622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1841715) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409055) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.7715001) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(-2.2489108) q[2];
sx q[2];
rz(-2.953139) q[2];
sx q[2];
rz(2.1072731) q[2];
rz(2.6876642) q[3];
sx q[3];
rz(-2.6700927) q[3];
sx q[3];
rz(2.3910458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];