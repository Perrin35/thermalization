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
rz(-0.79688537) q[0];
sx q[0];
rz(-1.4112043) q[0];
sx q[0];
rz(-2.222173) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(3.0923831) q[1];
sx q[1];
rz(10.13848) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7763166) q[0];
sx q[0];
rz(-2.8284024) q[0];
sx q[0];
rz(-1.0745144) q[0];
rz(-pi) q[1];
rz(3.0786985) q[2];
sx q[2];
rz(-1.4425292) q[2];
sx q[2];
rz(-0.72606444) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7138243) q[1];
sx q[1];
rz(-1.5694027) q[1];
sx q[1];
rz(-0.45098253) q[1];
rz(-1.6380129) q[3];
sx q[3];
rz(-0.47223642) q[3];
sx q[3];
rz(-0.41154644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.087622341) q[2];
sx q[2];
rz(-0.38603187) q[2];
sx q[2];
rz(2.7719882) q[2];
rz(1.1450279) q[3];
sx q[3];
rz(-1.4729045) q[3];
sx q[3];
rz(2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358148) q[0];
sx q[0];
rz(-1.4666297) q[0];
sx q[0];
rz(-0.57304397) q[0];
rz(2.0447958) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(-0.47164741) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1435463) q[0];
sx q[0];
rz(-0.61206383) q[0];
sx q[0];
rz(1.7630811) q[0];
x q[1];
rz(-0.71895968) q[2];
sx q[2];
rz(-2.3078354) q[2];
sx q[2];
rz(-1.3477088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1106893) q[1];
sx q[1];
rz(-0.86319268) q[1];
sx q[1];
rz(-2.3861107) q[1];
rz(0.40741323) q[3];
sx q[3];
rz(-1.0828567) q[3];
sx q[3];
rz(0.98516243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7219217) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(1.088885) q[2];
rz(1.0645083) q[3];
sx q[3];
rz(-2.0239315) q[3];
sx q[3];
rz(-2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0629145) q[0];
sx q[0];
rz(-2.9525472) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(0.71006376) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(0.56627083) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44646548) q[0];
sx q[0];
rz(-2.4003865) q[0];
sx q[0];
rz(-1.5542073) q[0];
x q[1];
rz(-2.5336419) q[2];
sx q[2];
rz(-2.3312116) q[2];
sx q[2];
rz(1.9160401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8108845) q[1];
sx q[1];
rz(-0.20527923) q[1];
sx q[1];
rz(0.42930557) q[1];
rz(-2.8763883) q[3];
sx q[3];
rz(-1.5062766) q[3];
sx q[3];
rz(1.8876739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.337734) q[2];
sx q[2];
rz(-0.89581076) q[2];
sx q[2];
rz(-3.1324978) q[2];
rz(-0.13633063) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(1.2135308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.424054) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(0.16615443) q[0];
rz(-2.0280139) q[1];
sx q[1];
rz(-0.49003092) q[1];
sx q[1];
rz(-2.9794433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067569392) q[0];
sx q[0];
rz(-2.9386006) q[0];
sx q[0];
rz(2.0503886) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5954564) q[2];
sx q[2];
rz(-2.1196105) q[2];
sx q[2];
rz(-1.2545214) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.051013324) q[1];
sx q[1];
rz(-0.83003488) q[1];
sx q[1];
rz(2.6464737) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1433853) q[3];
sx q[3];
rz(-0.95164883) q[3];
sx q[3];
rz(2.038508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1589511) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(-2.6137433) q[2];
rz(-0.85150254) q[3];
sx q[3];
rz(-2.8179171) q[3];
sx q[3];
rz(0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714137) q[0];
sx q[0];
rz(-1.5988007) q[0];
sx q[0];
rz(-0.80108368) q[0];
rz(-0.78394765) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(0.98091006) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.635467) q[0];
sx q[0];
rz(-2.6775595) q[0];
sx q[0];
rz(-2.6543174) q[0];
rz(-pi) q[1];
rz(1.0641685) q[2];
sx q[2];
rz(-1.9031823) q[2];
sx q[2];
rz(1.4719964) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0438761) q[1];
sx q[1];
rz(-1.4773158) q[1];
sx q[1];
rz(3.082117) q[1];
x q[2];
rz(-0.24215908) q[3];
sx q[3];
rz(-2.1091614) q[3];
sx q[3];
rz(0.70555609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.058978) q[2];
sx q[2];
rz(-1.5951944) q[2];
sx q[2];
rz(-0.8832461) q[2];
rz(2.9412681) q[3];
sx q[3];
rz(-0.89503461) q[3];
sx q[3];
rz(1.0909572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9969295) q[0];
sx q[0];
rz(-2.0893593) q[0];
sx q[0];
rz(-2.1494179) q[0];
rz(0.80157533) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(1.7209524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376037) q[0];
sx q[0];
rz(-2.1891356) q[0];
sx q[0];
rz(1.9644587) q[0];
x q[1];
rz(2.1031441) q[2];
sx q[2];
rz(-1.4803998) q[2];
sx q[2];
rz(2.9351075) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18773026) q[1];
sx q[1];
rz(-1.8988601) q[1];
sx q[1];
rz(-0.66098722) q[1];
rz(-0.037461683) q[3];
sx q[3];
rz(-1.6679296) q[3];
sx q[3];
rz(-2.284632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27986032) q[2];
sx q[2];
rz(-1.7143152) q[2];
sx q[2];
rz(-2.6045065) q[2];
rz(-3.1308657) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(-0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.895973) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(0.33988345) q[0];
rz(1.8396359) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(-1.9702912) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9635668) q[0];
sx q[0];
rz(-1.2739355) q[0];
sx q[0];
rz(1.6087449) q[0];
rz(2.3717959) q[2];
sx q[2];
rz(-0.30210051) q[2];
sx q[2];
rz(-0.090902791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.724546) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(-1.538365) q[1];
x q[2];
rz(0.62452353) q[3];
sx q[3];
rz(-1.4281264) q[3];
sx q[3];
rz(2.4017911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12477144) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(-0.67145124) q[2];
rz(2.6050383) q[3];
sx q[3];
rz(-2.1814929) q[3];
sx q[3];
rz(1.051739) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272098) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(0.89624727) q[0];
rz(-2.8889636) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(-2.0910697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80879656) q[0];
sx q[0];
rz(-2.0967297) q[0];
sx q[0];
rz(1.512708) q[0];
rz(-2.4765313) q[2];
sx q[2];
rz(-1.4675234) q[2];
sx q[2];
rz(-0.85368431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0922565) q[1];
sx q[1];
rz(-2.6529556) q[1];
sx q[1];
rz(-1.8948003) q[1];
x q[2];
rz(-3.0862634) q[3];
sx q[3];
rz(-2.3829975) q[3];
sx q[3];
rz(0.38842312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5184021) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(1.5790348) q[2];
rz(1.3671499) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(-0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52687454) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(0.38994625) q[0];
rz(-0.82603106) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(-1.0521851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5717333) q[0];
sx q[0];
rz(-0.89825199) q[0];
sx q[0];
rz(2.74617) q[0];
rz(-0.39689831) q[2];
sx q[2];
rz(-1.4087311) q[2];
sx q[2];
rz(-1.6610749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4349367) q[1];
sx q[1];
rz(-1.2806007) q[1];
sx q[1];
rz(1.1654084) q[1];
x q[2];
rz(-2.6652419) q[3];
sx q[3];
rz(-2.1423376) q[3];
sx q[3];
rz(-1.0494572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(-2.6123987) q[2];
rz(-0.75096327) q[3];
sx q[3];
rz(-0.96143985) q[3];
sx q[3];
rz(-2.9474337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26218364) q[0];
sx q[0];
rz(-0.59469596) q[0];
sx q[0];
rz(1.9770812) q[0];
rz(1.8544082) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(2.6374292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39162779) q[0];
sx q[0];
rz(-2.7996783) q[0];
sx q[0];
rz(0.86655273) q[0];
rz(-pi) q[1];
x q[1];
rz(1.590056) q[2];
sx q[2];
rz(-1.9474721) q[2];
sx q[2];
rz(0.45805675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6109687) q[1];
sx q[1];
rz(-1.8357616) q[1];
sx q[1];
rz(0.19719657) q[1];
rz(1.3066585) q[3];
sx q[3];
rz(-2.4396601) q[3];
sx q[3];
rz(0.86938996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21504271) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(-3.1070993) q[2];
rz(1.6132332) q[3];
sx q[3];
rz(-2.4210763) q[3];
sx q[3];
rz(0.31149402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(2.9813949) q[0];
sx q[0];
rz(-1.5934516) q[0];
sx q[0];
rz(-0.39504575) q[0];
rz(2.1389217) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(-0.38242292) q[2];
sx q[2];
rz(-2.6527846) q[2];
sx q[2];
rz(-2.6058886) q[2];
rz(-0.32714365) q[3];
sx q[3];
rz(-1.6302623) q[3];
sx q[3];
rz(0.74360256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
