OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(-2.350783) q[0];
sx q[0];
rz(2.8074582) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89361184) q[0];
sx q[0];
rz(-0.75151822) q[0];
sx q[0];
rz(0.052608842) q[0];
rz(-2.676744) q[2];
sx q[2];
rz(-1.5343622) q[2];
sx q[2];
rz(0.51793098) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7114746) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-0.75573604) q[1];
rz(-1.7928042) q[3];
sx q[3];
rz(-1.3862002) q[3];
sx q[3];
rz(2.8649462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-2.4770588) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537271) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4336006) q[0];
sx q[0];
rz(-3.1118244) q[0];
sx q[0];
rz(-3.004651) q[0];
rz(-2.5820929) q[2];
sx q[2];
rz(-1.3845978) q[2];
sx q[2];
rz(2.6452243) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2991997) q[1];
sx q[1];
rz(-2.0626915) q[1];
sx q[1];
rz(2.030034) q[1];
rz(-pi) q[2];
rz(1.6658695) q[3];
sx q[3];
rz(-0.94082309) q[3];
sx q[3];
rz(2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7188321) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(-2.9348532) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3690255) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(0.47779045) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-2.7405222) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7558407) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(1.8833478) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7602073) q[2];
sx q[2];
rz(-0.93511287) q[2];
sx q[2];
rz(2.1798101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9123147) q[1];
sx q[1];
rz(-0.73017987) q[1];
sx q[1];
rz(-0.0095403949) q[1];
rz(-pi) q[2];
rz(-2.0173666) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(-2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(0.25340432) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93638203) q[0];
sx q[0];
rz(-1.0316327) q[0];
sx q[0];
rz(-3.0217231) q[0];
rz(-pi) q[1];
x q[1];
rz(0.01732145) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(0.7893562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.036382347) q[1];
sx q[1];
rz(-1.6670334) q[1];
sx q[1];
rz(2.5114245) q[1];
rz(3.0074189) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(-1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0306586) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(0.91286719) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(1.0095989) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(0.062967904) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-2.6834992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040745) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(-3.1361561) q[0];
x q[1];
rz(2.6890254) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(2.3134311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3921567) q[1];
sx q[1];
rz(-1.6877618) q[1];
sx q[1];
rz(2.3163296) q[1];
rz(-pi) q[2];
rz(2.0714508) q[3];
sx q[3];
rz(-1.6541012) q[3];
sx q[3];
rz(-1.7617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5106817) q[0];
sx q[0];
rz(-1.1817389) q[0];
sx q[0];
rz(-0.77581497) q[0];
x q[1];
rz(0.74395545) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(0.26847408) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8900745) q[1];
sx q[1];
rz(-1.8971844) q[1];
sx q[1];
rz(0.20784394) q[1];
x q[2];
rz(-1.1250161) q[3];
sx q[3];
rz(-2.6187754) q[3];
sx q[3];
rz(1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1798114) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(2.069058) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(2.0297208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5892964) q[0];
sx q[0];
rz(-2.6483905) q[0];
sx q[0];
rz(-1.769355) q[0];
x q[1];
rz(-0.53084897) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(1.2209148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48549451) q[1];
sx q[1];
rz(-1.494325) q[1];
sx q[1];
rz(-0.2904201) q[1];
rz(-pi) q[2];
rz(0.38970077) q[3];
sx q[3];
rz(-2.0412363) q[3];
sx q[3];
rz(-2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(1.1676577) q[2];
rz(-1.6052823) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(1.4060098) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-2.3866167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687815) q[0];
sx q[0];
rz(-0.73505721) q[0];
sx q[0];
rz(1.789577) q[0];
x q[1];
rz(-0.72889363) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(0.58065562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4432115) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(-2.5618782) q[1];
rz(-pi) q[2];
rz(-2.2488392) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(1.0196109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76453152) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.586097) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(-3.1220904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855904) q[0];
sx q[0];
rz(-1.7770355) q[0];
sx q[0];
rz(0.081736728) q[0];
rz(1.3703913) q[2];
sx q[2];
rz(-1.7192171) q[2];
sx q[2];
rz(-2.5401149) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86917415) q[1];
sx q[1];
rz(-0.65009102) q[1];
sx q[1];
rz(-1.8723349) q[1];
x q[2];
rz(-2.9762514) q[3];
sx q[3];
rz(-1.2623566) q[3];
sx q[3];
rz(-0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(-1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(-1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.9932995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023037993) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(-1.5492357) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7933153) q[2];
sx q[2];
rz(-0.95138022) q[2];
sx q[2];
rz(-2.7811108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6069968) q[1];
sx q[1];
rz(-1.7417223) q[1];
sx q[1];
rz(-2.2439438) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(-1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(-3.0129516) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-2.1784492) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(2.6504436) q[2];
sx q[2];
rz(-0.72577234) q[2];
sx q[2];
rz(2.5189248) q[2];
rz(-2.1309489) q[3];
sx q[3];
rz(-1.5394566) q[3];
sx q[3];
rz(1.4395366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
