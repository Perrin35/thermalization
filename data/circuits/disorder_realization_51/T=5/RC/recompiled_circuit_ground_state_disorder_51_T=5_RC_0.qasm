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
rz(-0.33101141) q[1];
sx q[1];
rz(1.9537227) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5802104) q[0];
sx q[0];
rz(-2.8299709) q[0];
sx q[0];
rz(0.28579692) q[0];
x q[1];
rz(-1.3224363) q[2];
sx q[2];
rz(-2.6018086) q[2];
sx q[2];
rz(-2.7960461) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8459699) q[1];
sx q[1];
rz(-0.69087183) q[1];
sx q[1];
rz(0.9102896) q[1];
rz(-2.5099947) q[3];
sx q[3];
rz(-1.609841) q[3];
sx q[3];
rz(-1.4698879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0821705) q[2];
sx q[2];
rz(-2.1990081) q[2];
sx q[2];
rz(-1.4284632) q[2];
rz(-1.3905455) q[3];
sx q[3];
rz(-2.3640552) q[3];
sx q[3];
rz(2.9318504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9222337) q[0];
sx q[0];
rz(-1.5809504) q[0];
sx q[0];
rz(-1.6844164) q[0];
rz(-0.68151418) q[1];
sx q[1];
rz(-2.4512955) q[1];
sx q[1];
rz(2.1843279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0326704) q[0];
sx q[0];
rz(-1.3650948) q[0];
sx q[0];
rz(0.47423307) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25485867) q[2];
sx q[2];
rz(-0.83410848) q[2];
sx q[2];
rz(0.88803776) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87220472) q[1];
sx q[1];
rz(-1.1630668) q[1];
sx q[1];
rz(3.0881626) q[1];
rz(1.6396585) q[3];
sx q[3];
rz(-1.5605697) q[3];
sx q[3];
rz(0.49730147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7822781) q[2];
sx q[2];
rz(-2.094153) q[2];
sx q[2];
rz(-0.90685833) q[2];
rz(2.4685229) q[3];
sx q[3];
rz(-2.7570717) q[3];
sx q[3];
rz(-1.2955906) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(-1.9799045) q[1];
sx q[1];
rz(-3.0414157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4456383) q[0];
sx q[0];
rz(-1.7036519) q[0];
sx q[0];
rz(-1.8376875) q[0];
x q[1];
rz(-0.64073784) q[2];
sx q[2];
rz(-2.5362439) q[2];
sx q[2];
rz(-1.6289161) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.924739) q[1];
sx q[1];
rz(-1.7983872) q[1];
sx q[1];
rz(0.74149577) q[1];
rz(-2.4002714) q[3];
sx q[3];
rz(-0.59494114) q[3];
sx q[3];
rz(-0.5812318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.748041) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(1.0432517) q[2];
rz(2.4539963) q[3];
sx q[3];
rz(-0.50191534) q[3];
sx q[3];
rz(2.8858394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2954751) q[0];
sx q[0];
rz(-0.44545528) q[0];
sx q[0];
rz(-1.0365781) q[0];
rz(1.2713894) q[1];
sx q[1];
rz(-0.82415736) q[1];
sx q[1];
rz(-1.8545256) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6223479) q[0];
sx q[0];
rz(-0.96327269) q[0];
sx q[0];
rz(0.95444478) q[0];
x q[1];
rz(2.4802568) q[2];
sx q[2];
rz(-1.2244389) q[2];
sx q[2];
rz(-2.6544184) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27582622) q[1];
sx q[1];
rz(-2.4600907) q[1];
sx q[1];
rz(2.0382512) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.014145) q[3];
sx q[3];
rz(-1.2928767) q[3];
sx q[3];
rz(1.2025361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6048772) q[2];
sx q[2];
rz(-2.776919) q[2];
sx q[2];
rz(3.1216915) q[2];
rz(2.5455635) q[3];
sx q[3];
rz(-0.28087956) q[3];
sx q[3];
rz(1.9635487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.48536456) q[0];
sx q[0];
rz(-1.8940268) q[0];
sx q[0];
rz(1.9670991) q[0];
rz(-2.8354722) q[1];
sx q[1];
rz(-2.0968292) q[1];
sx q[1];
rz(-1.7721446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1533617) q[0];
sx q[0];
rz(-1.5772468) q[0];
sx q[0];
rz(3.1383443) q[0];
x q[1];
rz(0.17307504) q[2];
sx q[2];
rz(-2.1622116) q[2];
sx q[2];
rz(-1.8449291) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7983129) q[1];
sx q[1];
rz(-1.7213744) q[1];
sx q[1];
rz(1.2231138) q[1];
rz(-1.1569275) q[3];
sx q[3];
rz(-1.0936495) q[3];
sx q[3];
rz(-0.26354313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29641637) q[2];
sx q[2];
rz(-0.74411074) q[2];
sx q[2];
rz(2.3386492) q[2];
rz(2.0261649) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(-1.7723134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(0.95336103) q[0];
sx q[0];
rz(-0.53934923) q[0];
sx q[0];
rz(3.0292335) q[0];
rz(2.864783) q[1];
sx q[1];
rz(-1.1667629) q[1];
sx q[1];
rz(-0.64754957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4741199) q[0];
sx q[0];
rz(-1.9074567) q[0];
sx q[0];
rz(0.99322999) q[0];
rz(-pi) q[1];
rz(-0.55988257) q[2];
sx q[2];
rz(-0.75046021) q[2];
sx q[2];
rz(0.66302887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12880691) q[1];
sx q[1];
rz(-1.4527078) q[1];
sx q[1];
rz(-0.4507555) q[1];
rz(-pi) q[2];
rz(-2.9660951) q[3];
sx q[3];
rz(-1.3385404) q[3];
sx q[3];
rz(-1.0009223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4623146) q[2];
sx q[2];
rz(-2.9913112) q[2];
sx q[2];
rz(-1.798299) q[2];
rz(0.51370931) q[3];
sx q[3];
rz(-1.4880344) q[3];
sx q[3];
rz(-0.59521365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3419471) q[0];
sx q[0];
rz(-1.8193614) q[0];
sx q[0];
rz(0.42022002) q[0];
rz(1.3601903) q[1];
sx q[1];
rz(-1.1667292) q[1];
sx q[1];
rz(-0.80418783) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.490896) q[0];
sx q[0];
rz(-2.4593818) q[0];
sx q[0];
rz(-2.2006486) q[0];
rz(-pi) q[1];
rz(1.2857513) q[2];
sx q[2];
rz(-2.3346296) q[2];
sx q[2];
rz(-2.3318044) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.907885) q[1];
sx q[1];
rz(-1.4236439) q[1];
sx q[1];
rz(0.99214913) q[1];
x q[2];
rz(-1.4057833) q[3];
sx q[3];
rz(-1.4632153) q[3];
sx q[3];
rz(0.089311205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5516025) q[2];
sx q[2];
rz(-2.1253822) q[2];
sx q[2];
rz(1.1691947) q[2];
rz(-2.9679838) q[3];
sx q[3];
rz(-2.2326525) q[3];
sx q[3];
rz(-1.199523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6731897) q[0];
sx q[0];
rz(-2.115374) q[0];
sx q[0];
rz(-1.9386559) q[0];
rz(2.9007593) q[1];
sx q[1];
rz(-2.2151561) q[1];
sx q[1];
rz(-2.208272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77887452) q[0];
sx q[0];
rz(-1.2347455) q[0];
sx q[0];
rz(1.1707992) q[0];
rz(1.9674932) q[2];
sx q[2];
rz(-0.60985127) q[2];
sx q[2];
rz(1.6624474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49737469) q[1];
sx q[1];
rz(-1.3709732) q[1];
sx q[1];
rz(-2.7544061) q[1];
x q[2];
rz(-0.15214129) q[3];
sx q[3];
rz(-2.1295432) q[3];
sx q[3];
rz(-1.5882176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11747083) q[2];
sx q[2];
rz(-1.4535707) q[2];
sx q[2];
rz(-3.0061099) q[2];
rz(0.99902117) q[3];
sx q[3];
rz(-2.8935367) q[3];
sx q[3];
rz(2.5854056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713772) q[0];
sx q[0];
rz(-2.0129634) q[0];
sx q[0];
rz(1.7673329) q[0];
rz(-1.0820092) q[1];
sx q[1];
rz(-0.78649414) q[1];
sx q[1];
rz(-1.5793922) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0572206) q[0];
sx q[0];
rz(-1.6037774) q[0];
sx q[0];
rz(2.1103561) q[0];
x q[1];
rz(-1.1861997) q[2];
sx q[2];
rz(-2.0973258) q[2];
sx q[2];
rz(-2.1411454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7074896) q[1];
sx q[1];
rz(-1.130778) q[1];
sx q[1];
rz(2.773519) q[1];
rz(-pi) q[2];
rz(-2.0252941) q[3];
sx q[3];
rz(-2.4387932) q[3];
sx q[3];
rz(-1.9925607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18206) q[2];
sx q[2];
rz(-1.5675194) q[2];
sx q[2];
rz(-0.62425295) q[2];
rz(-0.67638451) q[3];
sx q[3];
rz(-1.1712149) q[3];
sx q[3];
rz(-2.306365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312209) q[0];
sx q[0];
rz(-1.2483163) q[0];
sx q[0];
rz(-1.7927908) q[0];
rz(-0.30385941) q[1];
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
rz(-2.9187885) q[0];
sx q[0];
rz(-2.3979146) q[0];
sx q[0];
rz(-0.2144465) q[0];
rz(-pi) q[1];
rz(-0.10909144) q[2];
sx q[2];
rz(-1.8938091) q[2];
sx q[2];
rz(-3.035088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59047952) q[1];
sx q[1];
rz(-2.5462228) q[1];
sx q[1];
rz(-0.84748419) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56352625) q[3];
sx q[3];
rz(-2.1057662) q[3];
sx q[3];
rz(-0.68171147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0554241) q[2];
sx q[2];
rz(-2.7358027) q[2];
sx q[2];
rz(2.1093192) q[2];
rz(-1.0885193) q[3];
sx q[3];
rz(-1.6531331) q[3];
sx q[3];
rz(-2.3672339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.1879723) q[0];
sx q[0];
rz(-0.80360501) q[0];
sx q[0];
rz(-3.0927717) q[0];
rz(1.8232518) q[1];
sx q[1];
rz(-1.7346458) q[1];
sx q[1];
rz(0.55191747) q[1];
rz(1.8404519) q[2];
sx q[2];
rz(-1.9741304) q[2];
sx q[2];
rz(-1.3000549) q[2];
rz(-2.0069176) q[3];
sx q[3];
rz(-1.5831309) q[3];
sx q[3];
rz(2.3203608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
