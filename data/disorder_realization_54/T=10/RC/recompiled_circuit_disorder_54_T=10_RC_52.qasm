OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78046507) q[0];
sx q[0];
rz(-2.4523003) q[0];
sx q[0];
rz(0.33049345) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5753484) q[0];
sx q[0];
rz(-0.69501221) q[0];
sx q[0];
rz(1.1027176) q[0];
rz(-pi) q[1];
rz(-1.66967) q[2];
sx q[2];
rz(-0.31497248) q[2];
sx q[2];
rz(1.0613943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82422772) q[1];
sx q[1];
rz(-1.0365651) q[1];
sx q[1];
rz(0.51170106) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0508918) q[3];
sx q[3];
rz(-1.2399779) q[3];
sx q[3];
rz(-1.0370805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4101397) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(1.6072134) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(2.5193135) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(2.2252749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.672357) q[0];
sx q[0];
rz(-1.0807481) q[0];
sx q[0];
rz(-2.310918) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4257405) q[2];
sx q[2];
rz(-2.3999891) q[2];
sx q[2];
rz(1.2765826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9333222) q[1];
sx q[1];
rz(-1.9083438) q[1];
sx q[1];
rz(-0.34668215) q[1];
rz(-0.022577062) q[3];
sx q[3];
rz(-1.1728247) q[3];
sx q[3];
rz(-0.28385362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24094412) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-2.3201578) q[2];
rz(-3.1243096) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(-1.0082555) q[0];
rz(-0.035765212) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(2.6170513) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4315955) q[0];
sx q[0];
rz(-1.734373) q[0];
sx q[0];
rz(-0.030199108) q[0];
rz(-pi) q[1];
rz(1.8345941) q[2];
sx q[2];
rz(-0.909415) q[2];
sx q[2];
rz(1.7922572) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14568612) q[1];
sx q[1];
rz(-1.5090764) q[1];
sx q[1];
rz(1.5476336) q[1];
rz(1.1509622) q[3];
sx q[3];
rz(-1.435278) q[3];
sx q[3];
rz(-0.47437048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77164578) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(-1.4952205) q[2];
rz(1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33048531) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-0.088949732) q[0];
rz(-2.6308909) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9289963) q[0];
sx q[0];
rz(-1.4082452) q[0];
sx q[0];
rz(-2.8430804) q[0];
rz(-2.5783522) q[2];
sx q[2];
rz(-13/(3*pi)) q[2];
sx q[2];
rz(-2.9134977) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62024414) q[1];
sx q[1];
rz(-2.2921831) q[1];
sx q[1];
rz(-1.1299302) q[1];
rz(-pi) q[2];
rz(-2.6570286) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(1.3467005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4758063) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(0.62292567) q[2];
rz(2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(-1.6437644) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4125815) q[0];
sx q[0];
rz(-2.3238365) q[0];
sx q[0];
rz(-1.9304995) q[0];
x q[1];
rz(-1.6696817) q[2];
sx q[2];
rz(-0.77557287) q[2];
sx q[2];
rz(-1.0727739) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0651107) q[1];
sx q[1];
rz(-1.5434192) q[1];
sx q[1];
rz(2.4379424) q[1];
rz(2.6260914) q[3];
sx q[3];
rz(-1.7121797) q[3];
sx q[3];
rz(-0.055671234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(1.5636469) q[2];
rz(2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6435796) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(-1.4354338) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5481789) q[0];
sx q[0];
rz(-0.81028623) q[0];
sx q[0];
rz(0.20472783) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4749182) q[2];
sx q[2];
rz(-2.5898858) q[2];
sx q[2];
rz(-0.72223896) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1286436) q[1];
sx q[1];
rz(-1.9777021) q[1];
sx q[1];
rz(-0.13601555) q[1];
rz(-pi) q[2];
rz(-3.1141838) q[3];
sx q[3];
rz(-1.8084744) q[3];
sx q[3];
rz(0.72565597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0662213) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(-0.071468778) q[2];
rz(-1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8554095) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(2.563971) q[0];
rz(1.2795992) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-2.1320027) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6972835) q[0];
sx q[0];
rz(-0.59519207) q[0];
sx q[0];
rz(-0.62023456) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7570653) q[2];
sx q[2];
rz(-1.782801) q[2];
sx q[2];
rz(0.42696135) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13751444) q[1];
sx q[1];
rz(-2.5234748) q[1];
sx q[1];
rz(2.532258) q[1];
rz(2.5261577) q[3];
sx q[3];
rz(-2.2781567) q[3];
sx q[3];
rz(0.59188852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(-0.47232929) q[3];
sx q[3];
rz(-1.6475369) q[3];
sx q[3];
rz(1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(2.902466) q[0];
rz(-0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-0.4447287) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0895773) q[0];
sx q[0];
rz(-1.9964295) q[0];
sx q[0];
rz(0.94469597) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93034805) q[2];
sx q[2];
rz(-2.1313416) q[2];
sx q[2];
rz(-0.7330187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4689444) q[1];
sx q[1];
rz(-2.4441507) q[1];
sx q[1];
rz(-1.3717321) q[1];
rz(-pi) q[2];
rz(2.8430311) q[3];
sx q[3];
rz(-0.42793722) q[3];
sx q[3];
rz(-0.096979389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(0.81531173) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(1.3172654) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(0.28717336) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-0.33219355) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591601) q[0];
sx q[0];
rz(-0.81572616) q[0];
sx q[0];
rz(2.0072323) q[0];
x q[1];
rz(-1.6261149) q[2];
sx q[2];
rz(-0.69960591) q[2];
sx q[2];
rz(-2.827022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5481944) q[1];
sx q[1];
rz(-1.8264923) q[1];
sx q[1];
rz(2.7138608) q[1];
rz(-2.513047) q[3];
sx q[3];
rz(-2.5124031) q[3];
sx q[3];
rz(-0.43720804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(-1.2906637) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-0.65657842) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(-1.0593876) q[0];
rz(-0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(-0.25451452) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6063639) q[0];
sx q[0];
rz(-1.6441206) q[0];
sx q[0];
rz(1.3396157) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8149257) q[2];
sx q[2];
rz(-2.0460528) q[2];
sx q[2];
rz(0.82820669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3855615) q[1];
sx q[1];
rz(-0.33547151) q[1];
sx q[1];
rz(-0.24753333) q[1];
x q[2];
rz(-2.9322846) q[3];
sx q[3];
rz(-0.44692398) q[3];
sx q[3];
rz(2.350654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.09482) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(2.0937031) q[2];
rz(1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027325252) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(2.7813773) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(-1.04503) q[2];
sx q[2];
rz(-1.3617931) q[2];
sx q[2];
rz(1.1522273) q[2];
rz(-0.23056728) q[3];
sx q[3];
rz(-2.1721526) q[3];
sx q[3];
rz(-0.9823907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
