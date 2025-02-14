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
rz(4.8719811) q[0];
sx q[0];
rz(13.48579) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(-0.04920955) q[1];
sx q[1];
rz(2.4278909) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27016923) q[0];
sx q[0];
rz(-1.7180301) q[0];
sx q[0];
rz(-1.8482313) q[0];
rz(-pi) q[1];
rz(3.0786985) q[2];
sx q[2];
rz(-1.6990635) q[2];
sx q[2];
rz(0.72606444) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42776838) q[1];
sx q[1];
rz(-1.57219) q[1];
sx q[1];
rz(0.45098253) q[1];
rz(-pi) q[2];
rz(-2.0421175) q[3];
sx q[3];
rz(-1.6013535) q[3];
sx q[3];
rz(1.9224642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0539703) q[2];
sx q[2];
rz(-0.38603187) q[2];
sx q[2];
rz(-0.36960441) q[2];
rz(-1.1450279) q[3];
sx q[3];
rz(-1.6686882) q[3];
sx q[3];
rz(-0.37231529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10577781) q[0];
sx q[0];
rz(-1.4666297) q[0];
sx q[0];
rz(-0.57304397) q[0];
rz(1.0967968) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(0.47164741) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7268611) q[0];
sx q[0];
rz(-1.4607753) q[0];
sx q[0];
rz(-2.1741406) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.198493) q[2];
sx q[2];
rz(-2.1617956) q[2];
sx q[2];
rz(-0.87795382) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0309033) q[1];
sx q[1];
rz(-2.2784) q[1];
sx q[1];
rz(-0.75548197) q[1];
rz(-pi) q[2];
rz(2.0949267) q[3];
sx q[3];
rz(-1.9283623) q[3];
sx q[3];
rz(2.3563354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4196709) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(-1.088885) q[2];
rz(-2.0770843) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(-0.3256807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0786781) q[0];
sx q[0];
rz(-2.9525472) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(-0.71006376) q[1];
sx q[1];
rz(-2.3080669) q[1];
sx q[1];
rz(-2.5753218) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46895263) q[0];
sx q[0];
rz(-0.8297161) q[0];
sx q[0];
rz(-3.1264105) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5336419) q[2];
sx q[2];
rz(-0.81038108) q[2];
sx q[2];
rz(-1.9160401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.10668392) q[1];
sx q[1];
rz(-1.3843754) q[1];
sx q[1];
rz(-1.4843462) q[1];
rz(-0.26520437) q[3];
sx q[3];
rz(-1.5062766) q[3];
sx q[3];
rz(1.2539188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8038586) q[2];
sx q[2];
rz(-0.89581076) q[2];
sx q[2];
rz(-3.1324978) q[2];
rz(0.13633063) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(1.9280619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.424054) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(-2.9754382) q[0];
rz(1.1135788) q[1];
sx q[1];
rz(-0.49003092) q[1];
sx q[1];
rz(0.16214935) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55570468) q[0];
sx q[0];
rz(-1.7506208) q[0];
sx q[0];
rz(3.0469045) q[0];
x q[1];
rz(-3.1012898) q[2];
sx q[2];
rz(-2.5922814) q[2];
sx q[2];
rz(-1.9343164) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0905793) q[1];
sx q[1];
rz(-2.3115578) q[1];
sx q[1];
rz(-0.49511893) q[1];
rz(-0.6643296) q[3];
sx q[3];
rz(-1.2264612) q[3];
sx q[3];
rz(0.20928247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98264155) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(0.52784935) q[2];
rz(2.2900901) q[3];
sx q[3];
rz(-2.8179171) q[3];
sx q[3];
rz(0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714137) q[0];
sx q[0];
rz(-1.5988007) q[0];
sx q[0];
rz(2.340509) q[0];
rz(-2.357645) q[1];
sx q[1];
rz(-2.3990217) q[1];
sx q[1];
rz(-2.1606826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5194395) q[0];
sx q[0];
rz(-1.781917) q[0];
sx q[0];
rz(-2.7252174) q[0];
rz(-pi) q[1];
rz(-0.95246117) q[2];
sx q[2];
rz(-2.5437069) q[2];
sx q[2];
rz(-0.43274227) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6740711) q[1];
sx q[1];
rz(-1.630012) q[1];
sx q[1];
rz(-1.4771512) q[1];
rz(-pi) q[2];
rz(-1.9526236) q[3];
sx q[3];
rz(-2.5562048) q[3];
sx q[3];
rz(-0.25661925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.058978) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(-0.8832461) q[2];
rz(-2.9412681) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(-2.0506355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14466318) q[0];
sx q[0];
rz(-2.0893593) q[0];
sx q[0];
rz(0.99217478) q[0];
rz(-2.3400173) q[1];
sx q[1];
rz(-1.8482607) q[1];
sx q[1];
rz(-1.7209524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4822749) q[0];
sx q[0];
rz(-2.4226696) q[0];
sx q[0];
rz(-2.6470967) q[0];
rz(-pi) q[1];
rz(-1.747521) q[2];
sx q[2];
rz(-2.6023539) q[2];
sx q[2];
rz(1.5163744) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1375422) q[1];
sx q[1];
rz(-0.95064771) q[1];
sx q[1];
rz(-1.9778848) q[1];
rz(-pi) q[2];
rz(-0.037461683) q[3];
sx q[3];
rz(-1.473663) q[3];
sx q[3];
rz(-0.85696062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8617323) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(0.53708616) q[2];
rz(-0.01072695) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(-0.33988345) q[0];
rz(-1.3019568) q[1];
sx q[1];
rz(-0.94568959) q[1];
sx q[1];
rz(1.9702912) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9635668) q[0];
sx q[0];
rz(-1.8676571) q[0];
sx q[0];
rz(1.5328477) q[0];
rz(2.9214462) q[2];
sx q[2];
rz(-1.7793806) q[2];
sx q[2];
rz(0.91503798) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4170467) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(1.6032277) q[1];
rz(-pi) q[2];
rz(1.3955437) q[3];
sx q[3];
rz(-0.95357663) q[3];
sx q[3];
rz(-2.4127236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.12477144) q[2];
sx q[2];
rz(-1.7243959) q[2];
sx q[2];
rz(0.67145124) q[2];
rz(2.6050383) q[3];
sx q[3];
rz(-2.1814929) q[3];
sx q[3];
rz(1.051739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8272098) q[0];
sx q[0];
rz(-2.0874513) q[0];
sx q[0];
rz(-2.2453454) q[0];
rz(0.25262901) q[1];
sx q[1];
rz(-2.2815506) q[1];
sx q[1];
rz(2.0910697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4481215) q[0];
sx q[0];
rz(-0.52883178) q[0];
sx q[0];
rz(0.099683925) q[0];
rz(-0.66506135) q[2];
sx q[2];
rz(-1.4675234) q[2];
sx q[2];
rz(-2.2879083) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0493361) q[1];
sx q[1];
rz(-0.4886371) q[1];
sx q[1];
rz(1.8948003) q[1];
rz(-pi) q[2];
rz(-0.055329247) q[3];
sx q[3];
rz(-2.3829975) q[3];
sx q[3];
rz(2.7531695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5184021) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(1.5625578) q[2];
rz(-1.3671499) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(-2.3794543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6147181) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(0.38994625) q[0];
rz(0.82603106) q[1];
sx q[1];
rz(-2.4470058) q[1];
sx q[1];
rz(2.0894076) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1601801) q[0];
sx q[0];
rz(-0.76428586) q[0];
sx q[0];
rz(2.0212964) q[0];
rz(-pi) q[1];
rz(0.39689831) q[2];
sx q[2];
rz(-1.4087311) q[2];
sx q[2];
rz(-1.4805178) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3996399) q[1];
sx q[1];
rz(-1.9583079) q[1];
sx q[1];
rz(-2.8273929) q[1];
rz(0.47635079) q[3];
sx q[3];
rz(-2.1423376) q[3];
sx q[3];
rz(-1.0494572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0803926) q[2];
sx q[2];
rz(-1.281851) q[2];
sx q[2];
rz(2.6123987) q[2];
rz(2.3906294) q[3];
sx q[3];
rz(-0.96143985) q[3];
sx q[3];
rz(-2.9474337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879409) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(-1.9770812) q[0];
rz(1.8544082) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(-0.50416344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7994294) q[0];
sx q[0];
rz(-1.3124046) q[0];
sx q[0];
rz(-2.9151205) q[0];
x q[1];
rz(0.37673925) q[2];
sx q[2];
rz(-1.5887056) q[2];
sx q[2];
rz(1.1056545) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2626942) q[1];
sx q[1];
rz(-2.8126908) q[1];
sx q[1];
rz(-0.94543381) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2553798) q[3];
sx q[3];
rz(-1.4014114) q[3];
sx q[3];
rz(0.90506314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9265499) q[2];
sx q[2];
rz(-1.5861009) q[2];
sx q[2];
rz(-3.1070993) q[2];
rz(-1.6132332) q[3];
sx q[3];
rz(-2.4210763) q[3];
sx q[3];
rz(2.8300986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9813949) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(1.002671) q[1];
sx q[1];
rz(-1.4613338) q[1];
sx q[1];
rz(2.7636539) q[1];
rz(-1.7667234) q[2];
sx q[2];
rz(-1.1200323) q[2];
sx q[2];
rz(-3.0333698) q[2];
rz(-1.6335842) q[3];
sx q[3];
rz(-1.2442524) q[3];
sx q[3];
rz(2.2942345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
