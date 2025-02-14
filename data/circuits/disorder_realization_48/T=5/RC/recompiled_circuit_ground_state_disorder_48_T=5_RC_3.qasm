OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(-0.25622955) q[0];
rz(0.37707075) q[1];
sx q[1];
rz(-1.7989676) q[1];
sx q[1];
rz(1.0599729) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7840756) q[0];
sx q[0];
rz(-1.4526443) q[0];
sx q[0];
rz(-1.4489001) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6057618) q[2];
sx q[2];
rz(-2.2630073) q[2];
sx q[2];
rz(-3.0193605) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89965883) q[1];
sx q[1];
rz(-2.8202758) q[1];
sx q[1];
rz(0.71996538) q[1];
rz(-pi) q[2];
rz(1.6698631) q[3];
sx q[3];
rz(-1.4270253) q[3];
sx q[3];
rz(-3.037775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18067351) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(1.4291576) q[2];
rz(2.5299431) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(3.070224) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83577689) q[0];
sx q[0];
rz(-0.80128765) q[0];
sx q[0];
rz(2.3961156) q[0];
rz(-1.9396793) q[1];
sx q[1];
rz(-1.8169836) q[1];
sx q[1];
rz(-1.036693) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7734402) q[0];
sx q[0];
rz(-2.7651909) q[0];
sx q[0];
rz(-1.8252116) q[0];
rz(-pi) q[1];
rz(0.027341893) q[2];
sx q[2];
rz(-0.76702709) q[2];
sx q[2];
rz(-2.1576144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0526317) q[1];
sx q[1];
rz(-0.76238576) q[1];
sx q[1];
rz(2.7247468) q[1];
x q[2];
rz(-0.58425957) q[3];
sx q[3];
rz(-1.7528894) q[3];
sx q[3];
rz(-2.3626773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0006973) q[2];
sx q[2];
rz(-1.8706198) q[2];
sx q[2];
rz(0.97266436) q[2];
rz(0.85727143) q[3];
sx q[3];
rz(-1.0664777) q[3];
sx q[3];
rz(-1.2681786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.641441) q[0];
sx q[0];
rz(-2.2730136) q[0];
sx q[0];
rz(-2.3386173) q[0];
rz(-2.4909486) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(0.21779901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9150216) q[0];
sx q[0];
rz(-0.75998291) q[0];
sx q[0];
rz(1.0081916) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62816633) q[2];
sx q[2];
rz(-1.3828613) q[2];
sx q[2];
rz(2.4357901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4068702) q[1];
sx q[1];
rz(-2.0592505) q[1];
sx q[1];
rz(-3.0367756) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5434615) q[3];
sx q[3];
rz(-2.7418828) q[3];
sx q[3];
rz(0.78228116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.190072) q[2];
sx q[2];
rz(-2.0900487) q[2];
sx q[2];
rz(1.5104843) q[2];
rz(-1.914628) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.889582) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(-0.66396436) q[0];
rz(-0.12531677) q[1];
sx q[1];
rz(-1.7071416) q[1];
sx q[1];
rz(-2.1024316) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17500278) q[0];
sx q[0];
rz(-0.010207264) q[0];
sx q[0];
rz(-0.27423476) q[0];
rz(-pi) q[1];
rz(-2.8739615) q[2];
sx q[2];
rz(-1.9534785) q[2];
sx q[2];
rz(0.022993372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81280983) q[1];
sx q[1];
rz(-2.610958) q[1];
sx q[1];
rz(-1.1483436) q[1];
x q[2];
rz(-1.6658721) q[3];
sx q[3];
rz(-2.0119442) q[3];
sx q[3];
rz(0.26218647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8404954) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(3.0636129) q[2];
rz(1.1178499) q[3];
sx q[3];
rz(-1.0406787) q[3];
sx q[3];
rz(2.1054721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2074821) q[0];
sx q[0];
rz(-0.20620646) q[0];
sx q[0];
rz(-0.41314405) q[0];
rz(-1.235599) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(-1.8280169) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13554561) q[0];
sx q[0];
rz(-2.1925547) q[0];
sx q[0];
rz(-0.22613392) q[0];
rz(0.35982015) q[2];
sx q[2];
rz(-1.1319379) q[2];
sx q[2];
rz(-0.0024589389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47258618) q[1];
sx q[1];
rz(-2.5372977) q[1];
sx q[1];
rz(-2.2766606) q[1];
rz(0.78727874) q[3];
sx q[3];
rz(-1.9036999) q[3];
sx q[3];
rz(-1.5008139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.95973394) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(2.6971297) q[2];
rz(1.9410761) q[3];
sx q[3];
rz(-1.8903172) q[3];
sx q[3];
rz(-3.0357231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7926517) q[0];
sx q[0];
rz(-0.3949202) q[0];
sx q[0];
rz(0.077202395) q[0];
rz(-0.96241799) q[1];
sx q[1];
rz(-1.9541157) q[1];
sx q[1];
rz(0.033800689) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12594906) q[0];
sx q[0];
rz(-1.9140049) q[0];
sx q[0];
rz(0.48114602) q[0];
x q[1];
rz(1.5977497) q[2];
sx q[2];
rz(-2.1183062) q[2];
sx q[2];
rz(1.6545878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9056184) q[1];
sx q[1];
rz(-1.8840569) q[1];
sx q[1];
rz(2.226023) q[1];
x q[2];
rz(0.77671364) q[3];
sx q[3];
rz(-1.4956105) q[3];
sx q[3];
rz(0.0089587072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9683711) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(0.25924337) q[2];
rz(2.1974468) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(1.3313782) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060870085) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(-0.61332214) q[0];
rz(-2.2382286) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(-2.9936252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67576236) q[0];
sx q[0];
rz(-2.0480541) q[0];
sx q[0];
rz(-0.81264074) q[0];
x q[1];
rz(1.5388707) q[2];
sx q[2];
rz(-1.1319926) q[2];
sx q[2];
rz(1.1915468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2206055) q[1];
sx q[1];
rz(-1.1488394) q[1];
sx q[1];
rz(0.76442952) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1359753) q[3];
sx q[3];
rz(-1.5907173) q[3];
sx q[3];
rz(-1.7651778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.938574) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(-0.99986783) q[2];
rz(-1.8262919) q[3];
sx q[3];
rz(-0.2838997) q[3];
sx q[3];
rz(0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1891747) q[0];
sx q[0];
rz(-0.84680951) q[0];
sx q[0];
rz(-2.1536105) q[0];
rz(-1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(0.16407897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9768499) q[0];
sx q[0];
rz(-1.688042) q[0];
sx q[0];
rz(1.4848723) q[0];
rz(-pi) q[1];
rz(2.0664178) q[2];
sx q[2];
rz(-1.3775423) q[2];
sx q[2];
rz(0.19710625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2365109) q[1];
sx q[1];
rz(-1.0781546) q[1];
sx q[1];
rz(0.0035159455) q[1];
x q[2];
rz(2.9470575) q[3];
sx q[3];
rz(-2.2119868) q[3];
sx q[3];
rz(-1.7879144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24171955) q[2];
sx q[2];
rz(-0.49979979) q[2];
sx q[2];
rz(2.0797753) q[2];
rz(-0.1344943) q[3];
sx q[3];
rz(-2.2420292) q[3];
sx q[3];
rz(0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6580842) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(0.2051951) q[0];
rz(-0.23458734) q[1];
sx q[1];
rz(-1.7154452) q[1];
sx q[1];
rz(-2.9451784) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72442833) q[0];
sx q[0];
rz(-2.1311893) q[0];
sx q[0];
rz(0.36720328) q[0];
rz(1.4142198) q[2];
sx q[2];
rz(-1.051826) q[2];
sx q[2];
rz(-1.165498) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3383693) q[1];
sx q[1];
rz(-1.3937989) q[1];
sx q[1];
rz(1.2425674) q[1];
rz(0.78077448) q[3];
sx q[3];
rz(-1.1765683) q[3];
sx q[3];
rz(-0.0782644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74497574) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(-1.1783925) q[2];
rz(-0.34934238) q[3];
sx q[3];
rz(-1.1290461) q[3];
sx q[3];
rz(2.0918625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0013393764) q[0];
sx q[0];
rz(-2.0297191) q[0];
sx q[0];
rz(-2.4358791) q[0];
rz(2.4440675) q[1];
sx q[1];
rz(-1.768528) q[1];
sx q[1];
rz(2.2360905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6827992) q[0];
sx q[0];
rz(-1.1371277) q[0];
sx q[0];
rz(-1.5894029) q[0];
rz(1.897325) q[2];
sx q[2];
rz(-1.9166319) q[2];
sx q[2];
rz(-2.2188733) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4982371) q[1];
sx q[1];
rz(-1.61803) q[1];
sx q[1];
rz(-0.31075732) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9344744) q[3];
sx q[3];
rz(-1.26059) q[3];
sx q[3];
rz(-2.9077173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2424348) q[2];
sx q[2];
rz(-2.0826715) q[2];
sx q[2];
rz(-2.4230797) q[2];
rz(2.4588623) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(3.0375286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8231507) q[0];
sx q[0];
rz(-2.5704076) q[0];
sx q[0];
rz(3.0388863) q[0];
rz(1.6406583) q[1];
sx q[1];
rz(-1.8713015) q[1];
sx q[1];
rz(-0.4458977) q[1];
rz(-0.6871625) q[2];
sx q[2];
rz(-2.9570815) q[2];
sx q[2];
rz(2.2312192) q[2];
rz(-2.4950776) q[3];
sx q[3];
rz(-0.91520354) q[3];
sx q[3];
rz(-1.147816) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
