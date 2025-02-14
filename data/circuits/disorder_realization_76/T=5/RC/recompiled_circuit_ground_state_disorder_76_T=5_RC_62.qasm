OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(-1.0816242) q[0];
sx q[0];
rz(2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(-1.7159599) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65042881) q[0];
sx q[0];
rz(-0.61507498) q[0];
sx q[0];
rz(-2.0838724) q[0];
rz(2.6447634) q[2];
sx q[2];
rz(-1.8112) q[2];
sx q[2];
rz(0.59076004) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5162075) q[1];
sx q[1];
rz(-2.974286) q[1];
sx q[1];
rz(1.0870666) q[1];
rz(-pi) q[2];
rz(3.0979068) q[3];
sx q[3];
rz(-2.2891392) q[3];
sx q[3];
rz(1.1506611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90929675) q[2];
sx q[2];
rz(-0.84315073) q[2];
sx q[2];
rz(-0.24093957) q[2];
rz(0.029189261) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(-1.1791112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772684) q[0];
sx q[0];
rz(-1.203953) q[0];
sx q[0];
rz(-0.48467317) q[0];
rz(-0.67972216) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(1.9814804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77185936) q[0];
sx q[0];
rz(-1.5892649) q[0];
sx q[0];
rz(-2.8375677) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76267879) q[2];
sx q[2];
rz(-0.87717036) q[2];
sx q[2];
rz(2.9528303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45534652) q[1];
sx q[1];
rz(-2.1619611) q[1];
sx q[1];
rz(1.2711729) q[1];
rz(-pi) q[2];
rz(-1.2413361) q[3];
sx q[3];
rz(-2.0700392) q[3];
sx q[3];
rz(1.3105621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.844187) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(2.6837132) q[2];
rz(2.0186021) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(-2.5751953) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26894012) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(3.1100682) q[0];
rz(-2.8541376) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.215975) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7087962) q[0];
sx q[0];
rz(-2.2026718) q[0];
sx q[0];
rz(-2.2461476) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38162614) q[2];
sx q[2];
rz(-0.51593057) q[2];
sx q[2];
rz(2.6037773) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79920125) q[1];
sx q[1];
rz(-1.2974129) q[1];
sx q[1];
rz(-2.8552516) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3787556) q[3];
sx q[3];
rz(-2.3194072) q[3];
sx q[3];
rz(1.1824106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.03881255) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(-2.3056324) q[2];
rz(-1.3577667) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(0.74762216) q[3];
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
rz(pi/2) q[0];
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
rz(-1.8589856) q[0];
sx q[0];
rz(-1.405412) q[0];
sx q[0];
rz(-0.080168515) q[0];
rz(-1.0824341) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(-1.8128043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4565312) q[0];
sx q[0];
rz(-2.079112) q[0];
sx q[0];
rz(-1.9446745) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8187739) q[2];
sx q[2];
rz(-1.7344966) q[2];
sx q[2];
rz(1.9893008) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.13212407) q[1];
sx q[1];
rz(-1.7156202) q[1];
sx q[1];
rz(2.0028466) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6065663) q[3];
sx q[3];
rz(-0.38540977) q[3];
sx q[3];
rz(0.42675323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7666011) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(0.90744606) q[2];
rz(0.67029101) q[3];
sx q[3];
rz(-0.82796103) q[3];
sx q[3];
rz(-1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9012673) q[0];
sx q[0];
rz(-0.87258744) q[0];
sx q[0];
rz(0.43148828) q[0];
rz(1.1314499) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(2.7010837) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1898855) q[0];
sx q[0];
rz(-1.110297) q[0];
sx q[0];
rz(3.0104464) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9856057) q[2];
sx q[2];
rz(-1.3813263) q[2];
sx q[2];
rz(1.0111601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8590282) q[1];
sx q[1];
rz(-2.8088114) q[1];
sx q[1];
rz(-3.0206693) q[1];
rz(-0.83310762) q[3];
sx q[3];
rz(-2.3596977) q[3];
sx q[3];
rz(1.9946919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12209192) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(0.41395536) q[2];
rz(1.7533938) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(3.0164914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269161) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(2.2077014) q[0];
rz(2.4712708) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(1.3353039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40437296) q[0];
sx q[0];
rz(-1.6707194) q[0];
sx q[0];
rz(-2.1735682) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15768965) q[2];
sx q[2];
rz(-1.5957498) q[2];
sx q[2];
rz(1.7316929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.851593) q[1];
sx q[1];
rz(-1.6834752) q[1];
sx q[1];
rz(-3.1188909) q[1];
x q[2];
rz(1.5621095) q[3];
sx q[3];
rz(-2.207805) q[3];
sx q[3];
rz(-0.52906936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2495217) q[2];
sx q[2];
rz(-0.27270174) q[2];
sx q[2];
rz(1.8708694) q[2];
rz(-1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(-0.52186596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524566) q[0];
sx q[0];
rz(-2.554775) q[0];
sx q[0];
rz(0.15368803) q[0];
rz(0.20248374) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(0.94857803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62371333) q[0];
sx q[0];
rz(-2.1800632) q[0];
sx q[0];
rz(1.3246956) q[0];
rz(-pi) q[1];
rz(-2.2099451) q[2];
sx q[2];
rz(-1.3077362) q[2];
sx q[2];
rz(-1.6073038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2413509) q[1];
sx q[1];
rz(-0.35120041) q[1];
sx q[1];
rz(2.3022149) q[1];
rz(-1.7465215) q[3];
sx q[3];
rz(-0.65726377) q[3];
sx q[3];
rz(-0.45096179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.004868) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(-0.0014121545) q[2];
rz(-0.062156113) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(2.1803161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75981265) q[0];
sx q[0];
rz(-2.3842922) q[0];
sx q[0];
rz(-1.9267474) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-2.6140723) q[1];
sx q[1];
rz(0.55799276) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7400949) q[0];
sx q[0];
rz(-1.8170274) q[0];
sx q[0];
rz(1.2241062) q[0];
rz(-pi) q[1];
rz(-2.9079014) q[2];
sx q[2];
rz(-0.56392852) q[2];
sx q[2];
rz(-0.83966161) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3000189) q[1];
sx q[1];
rz(-1.2928559) q[1];
sx q[1];
rz(-2.071408) q[1];
rz(-1.4766944) q[3];
sx q[3];
rz(-2.8633139) q[3];
sx q[3];
rz(-2.813193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5672292) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(-0.26926678) q[2];
rz(1.1632129) q[3];
sx q[3];
rz(-0.60267699) q[3];
sx q[3];
rz(-0.24063024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6729386) q[0];
sx q[0];
rz(-2.1042295) q[0];
sx q[0];
rz(1.0733676) q[0];
rz(-2.0138373) q[1];
sx q[1];
rz(-0.22722166) q[1];
sx q[1];
rz(-2.5849297) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4024324) q[0];
sx q[0];
rz(-2.7910829) q[0];
sx q[0];
rz(1.1465766) q[0];
rz(-pi) q[1];
rz(-0.9999335) q[2];
sx q[2];
rz(-2.0682671) q[2];
sx q[2];
rz(1.2116878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3240668) q[1];
sx q[1];
rz(-0.81302596) q[1];
sx q[1];
rz(-0.50285411) q[1];
x q[2];
rz(-1.6708371) q[3];
sx q[3];
rz(-1.3373858) q[3];
sx q[3];
rz(1.4518713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2399981) q[2];
sx q[2];
rz(-1.0866714) q[2];
sx q[2];
rz(1.979801) q[2];
rz(0.53317541) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(-2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.48851442) q[0];
sx q[0];
rz(-0.97452679) q[0];
sx q[0];
rz(-2.5467806) q[0];
rz(-2.5626903) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(-2.4822809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4730277) q[0];
sx q[0];
rz(-0.94946948) q[0];
sx q[0];
rz(-1.7354119) q[0];
x q[1];
rz(-2.5393007) q[2];
sx q[2];
rz(-2.4458439) q[2];
sx q[2];
rz(-0.93590036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6281471) q[1];
sx q[1];
rz(-1.7127348) q[1];
sx q[1];
rz(-1.3238841) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3385184) q[3];
sx q[3];
rz(-2.7702906) q[3];
sx q[3];
rz(0.40483958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79002964) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(-0.053038049) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(3.067335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2863083) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(0.098943624) q[1];
sx q[1];
rz(-1.3700486) q[1];
sx q[1];
rz(0.055421967) q[1];
rz(-2.7361353) q[2];
sx q[2];
rz(-1.5628855) q[2];
sx q[2];
rz(0.93456675) q[2];
rz(0.17114279) q[3];
sx q[3];
rz(-2.3928693) q[3];
sx q[3];
rz(0.00015043845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
