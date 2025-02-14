OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(2.0164665) q[0];
rz(-2.1195124) q[1];
sx q[1];
rz(-2.4740969) q[1];
sx q[1];
rz(-0.8134841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33685499) q[0];
sx q[0];
rz(-1.2053524) q[0];
sx q[0];
rz(2.9904537) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6935319) q[2];
sx q[2];
rz(-1.8984814) q[2];
sx q[2];
rz(1.1111914) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5578545) q[1];
sx q[1];
rz(-2.0101133) q[1];
sx q[1];
rz(-1.735461) q[1];
x q[2];
rz(-1.7343821) q[3];
sx q[3];
rz(-1.4814875) q[3];
sx q[3];
rz(-2.341604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(-0.92864621) q[2];
rz(-1.5993902) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9376675) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(0.12705886) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.3652722) q[1];
sx q[1];
rz(-0.7712706) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.49202) q[0];
sx q[0];
rz(-0.82081534) q[0];
sx q[0];
rz(1.0787021) q[0];
rz(-3.0808582) q[2];
sx q[2];
rz(-1.2337451) q[2];
sx q[2];
rz(3.0063546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88395547) q[1];
sx q[1];
rz(-1.9169382) q[1];
sx q[1];
rz(0.37948541) q[1];
rz(-pi) q[2];
rz(-2.1089606) q[3];
sx q[3];
rz(-2.0797044) q[3];
sx q[3];
rz(2.5050688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1473006) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(3.1331983) q[2];
rz(0.66347915) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10107772) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(-2.7000632) q[0];
rz(-2.1532374) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(3.006014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5268742) q[0];
sx q[0];
rz(-1.5784987) q[0];
sx q[0];
rz(0.016250261) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1889815) q[2];
sx q[2];
rz(-2.2509607) q[2];
sx q[2];
rz(1.176468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49040321) q[1];
sx q[1];
rz(-2.5156543) q[1];
sx q[1];
rz(2.2041026) q[1];
rz(-1.3531923) q[3];
sx q[3];
rz(-2.3953468) q[3];
sx q[3];
rz(-1.3775795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2077937) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(0.03820339) q[2];
rz(-2.6162052) q[3];
sx q[3];
rz(-2.5140258) q[3];
sx q[3];
rz(-0.6853404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66048375) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(2.5307181) q[0];
rz(-1.8065709) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(0.23922051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6253875) q[0];
sx q[0];
rz(-0.62407485) q[0];
sx q[0];
rz(2.027987) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47232136) q[2];
sx q[2];
rz(-1.2178253) q[2];
sx q[2];
rz(-1.8582839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4799616) q[1];
sx q[1];
rz(-1.1203655) q[1];
sx q[1];
rz(1.2110787) q[1];
rz(1.3705105) q[3];
sx q[3];
rz(-2.0242175) q[3];
sx q[3];
rz(-3.0016921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7188344) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(0.0017496721) q[2];
rz(-2.8478029) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(0.26688117) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9002429) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(-2.2985261) q[0];
rz(0.33755606) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(-1.5379803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6809236) q[0];
sx q[0];
rz(-2.3516555) q[0];
sx q[0];
rz(1.3000751) q[0];
x q[1];
rz(1.5174687) q[2];
sx q[2];
rz(-1.7449208) q[2];
sx q[2];
rz(-1.8873896) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5953491) q[1];
sx q[1];
rz(-1.8270632) q[1];
sx q[1];
rz(-1.3687737) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20974737) q[3];
sx q[3];
rz(-1.4919859) q[3];
sx q[3];
rz(2.1926853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43188492) q[2];
sx q[2];
rz(-2.0543435) q[2];
sx q[2];
rz(2.0951927) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-2.4304978) q[3];
sx q[3];
rz(-0.54767245) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(2.4936254) q[0];
rz(-1.3994392) q[1];
sx q[1];
rz(-1.6439227) q[1];
sx q[1];
rz(-1.8125777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8662939) q[0];
sx q[0];
rz(-1.4614551) q[0];
sx q[0];
rz(1.3491531) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7850661) q[2];
sx q[2];
rz(-1.2870711) q[2];
sx q[2];
rz(1.2905215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15550286) q[1];
sx q[1];
rz(-2.2515319) q[1];
sx q[1];
rz(2.7377241) q[1];
x q[2];
rz(-1.2043318) q[3];
sx q[3];
rz(-2.0559337) q[3];
sx q[3];
rz(-1.7462891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38137388) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(-1.5927429) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.2589688) q[3];
sx q[3];
rz(-0.86863345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233362) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(0.94183952) q[0];
rz(-2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-0.19217415) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1051467) q[0];
sx q[0];
rz(-0.22686401) q[0];
sx q[0];
rz(1.731621) q[0];
rz(1.6404387) q[2];
sx q[2];
rz(-0.52280871) q[2];
sx q[2];
rz(-1.4774587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3698764) q[1];
sx q[1];
rz(-1.2315005) q[1];
sx q[1];
rz(-2.3849065) q[1];
rz(-pi) q[2];
rz(-0.022054733) q[3];
sx q[3];
rz(-2.0682749) q[3];
sx q[3];
rz(-2.9201404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3840702) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(2.7080217) q[2];
rz(1.5844257) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(0.43237329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.3766377) q[0];
sx q[0];
rz(-1.7973416) q[0];
rz(-1.2035707) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(1.3628091) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840795) q[0];
sx q[0];
rz(-1.5410893) q[0];
sx q[0];
rz(-3.1196575) q[0];
rz(-pi) q[1];
rz(1.3819456) q[2];
sx q[2];
rz(-2.3879693) q[2];
sx q[2];
rz(0.058503956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4688022) q[1];
sx q[1];
rz(-1.6957449) q[1];
sx q[1];
rz(-2.9114257) q[1];
rz(-pi) q[2];
rz(2.2572956) q[3];
sx q[3];
rz(-1.3217759) q[3];
sx q[3];
rz(-0.62307138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5726996) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(0.9838689) q[2];
rz(-0.24108663) q[3];
sx q[3];
rz(-2.2086996) q[3];
sx q[3];
rz(-0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47138658) q[0];
sx q[0];
rz(-0.79600483) q[0];
sx q[0];
rz(2.8669226) q[0];
rz(-1.341691) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(2.0955657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957841) q[0];
sx q[0];
rz(-1.8950476) q[0];
sx q[0];
rz(-2.5155873) q[0];
x q[1];
rz(0.080950254) q[2];
sx q[2];
rz(-1.3737203) q[2];
sx q[2];
rz(2.112889) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.14635) q[1];
sx q[1];
rz(-2.1533794) q[1];
sx q[1];
rz(-1.6487153) q[1];
x q[2];
rz(2.0738129) q[3];
sx q[3];
rz(-1.5340337) q[3];
sx q[3];
rz(-3.0779148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.7509165) q[2];
sx q[2];
rz(0.39946237) q[2];
rz(-0.56600371) q[3];
sx q[3];
rz(-0.96650201) q[3];
sx q[3];
rz(2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-0.28867662) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(-0.5823108) q[0];
rz(-2.9583926) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(0.80642548) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9574653) q[0];
sx q[0];
rz(-0.13617198) q[0];
sx q[0];
rz(-2.2720112) q[0];
rz(-pi) q[1];
rz(2.8060032) q[2];
sx q[2];
rz(-1.8238471) q[2];
sx q[2];
rz(2.0584681) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6567071) q[1];
sx q[1];
rz(-2.2995512) q[1];
sx q[1];
rz(-1.2970112) q[1];
rz(-pi) q[2];
x q[2];
rz(1.459645) q[3];
sx q[3];
rz(-1.9160929) q[3];
sx q[3];
rz(-2.2636641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9091984) q[2];
sx q[2];
rz(-2.4269203) q[2];
sx q[2];
rz(2.3363028) q[2];
rz(-2.9946839) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-1.2522226) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(-1.5421142) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(0.11453828) q[2];
sx q[2];
rz(-2.2907612) q[2];
sx q[2];
rz(-1.9775122) q[2];
rz(2.6977758) q[3];
sx q[3];
rz(-2.4469821) q[3];
sx q[3];
rz(-1.2367005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
