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
rz(0.4515689) q[0];
sx q[0];
rz(-1.7234001) q[0];
sx q[0];
rz(-2.7127142) q[0];
rz(-0.14313993) q[1];
sx q[1];
rz(-2.0506471) q[1];
sx q[1];
rz(0.080088869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1983436) q[0];
sx q[0];
rz(-0.08028537) q[0];
sx q[0];
rz(-0.16585089) q[0];
rz(-1.3037138) q[2];
sx q[2];
rz(-0.3729698) q[2];
sx q[2];
rz(-2.3139062) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.061191843) q[1];
sx q[1];
rz(-2.3843003) q[1];
sx q[1];
rz(1.9593092) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4065885) q[3];
sx q[3];
rz(-0.99942151) q[3];
sx q[3];
rz(-0.77916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0668138) q[2];
sx q[2];
rz(-1.0513693) q[2];
sx q[2];
rz(-1.7195513) q[2];
rz(1.7285796) q[3];
sx q[3];
rz(-1.6001817) q[3];
sx q[3];
rz(-1.6933256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1220575) q[0];
sx q[0];
rz(-1.7449361) q[0];
sx q[0];
rz(-0.33357093) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-2.797778) q[1];
sx q[1];
rz(0.75210345) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010875794) q[0];
sx q[0];
rz(-2.1701522) q[0];
sx q[0];
rz(-2.5917971) q[0];
rz(1.5444438) q[2];
sx q[2];
rz(-1.7894723) q[2];
sx q[2];
rz(1.5062576) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.32622329) q[1];
sx q[1];
rz(-2.0953676) q[1];
sx q[1];
rz(-3.0154445) q[1];
rz(-1.7957572) q[3];
sx q[3];
rz(-1.3153425) q[3];
sx q[3];
rz(0.47522241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0636474) q[2];
sx q[2];
rz(-0.91219488) q[2];
sx q[2];
rz(0.65703854) q[2];
rz(-2.4296203) q[3];
sx q[3];
rz(-2.4896121) q[3];
sx q[3];
rz(0.74487346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4350568) q[0];
sx q[0];
rz(-0.61921316) q[0];
sx q[0];
rz(2.1001429) q[0];
rz(-0.374818) q[1];
sx q[1];
rz(-1.5033787) q[1];
sx q[1];
rz(-0.55694881) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0335046) q[0];
sx q[0];
rz(-1.9226941) q[0];
sx q[0];
rz(-0.75893388) q[0];
rz(-pi) q[1];
rz(1.7250437) q[2];
sx q[2];
rz(-0.45636141) q[2];
sx q[2];
rz(-0.92044324) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8136602) q[1];
sx q[1];
rz(-2.2371462) q[1];
sx q[1];
rz(-1.440669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3339504) q[3];
sx q[3];
rz(-2.2082553) q[3];
sx q[3];
rz(1.4887631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53691429) q[2];
sx q[2];
rz(-2.2173209) q[2];
sx q[2];
rz(0.43517819) q[2];
rz(-2.6168881) q[3];
sx q[3];
rz(-2.8080431) q[3];
sx q[3];
rz(-0.13478002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0404469) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(1.5144298) q[0];
rz(1.0487522) q[1];
sx q[1];
rz(-2.2188413) q[1];
sx q[1];
rz(1.705816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41489515) q[0];
sx q[0];
rz(-1.5279307) q[0];
sx q[0];
rz(2.9662762) q[0];
rz(2.5982791) q[2];
sx q[2];
rz(-1.8738973) q[2];
sx q[2];
rz(1.4599279) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2843316) q[1];
sx q[1];
rz(-1.5067717) q[1];
sx q[1];
rz(1.048719) q[1];
rz(-pi) q[2];
rz(0.34603186) q[3];
sx q[3];
rz(-1.1308987) q[3];
sx q[3];
rz(-1.7302681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40020308) q[2];
sx q[2];
rz(-1.783739) q[2];
sx q[2];
rz(0.025156585) q[2];
rz(1.4870421) q[3];
sx q[3];
rz(-2.7436723) q[3];
sx q[3];
rz(2.3253843) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095489278) q[0];
sx q[0];
rz(-0.64952055) q[0];
sx q[0];
rz(2.982614) q[0];
rz(-2.036463) q[1];
sx q[1];
rz(-2.4159894) q[1];
sx q[1];
rz(0.22430688) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0889322) q[0];
sx q[0];
rz(-1.972368) q[0];
sx q[0];
rz(-0.60202718) q[0];
rz(1.0105466) q[2];
sx q[2];
rz(-1.59861) q[2];
sx q[2];
rz(1.5956439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4233647) q[1];
sx q[1];
rz(-1.0486239) q[1];
sx q[1];
rz(-3.0084064) q[1];
x q[2];
rz(-2.9486233) q[3];
sx q[3];
rz(-2.5558839) q[3];
sx q[3];
rz(-1.0568413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9362681) q[2];
sx q[2];
rz(-1.779413) q[2];
sx q[2];
rz(-2.3929907) q[2];
rz(0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(-0.41142472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9801789) q[0];
sx q[0];
rz(-2.7730589) q[0];
sx q[0];
rz(-0.73032105) q[0];
rz(1.5018564) q[1];
sx q[1];
rz(-1.390099) q[1];
sx q[1];
rz(-2.9072442) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6662009) q[0];
sx q[0];
rz(-2.0115543) q[0];
sx q[0];
rz(-0.73385629) q[0];
rz(0.5989093) q[2];
sx q[2];
rz(-1.0157758) q[2];
sx q[2];
rz(2.2428494) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4428574) q[1];
sx q[1];
rz(-2.9720829) q[1];
sx q[1];
rz(-0.28556602) q[1];
rz(-pi) q[2];
rz(-1.1612915) q[3];
sx q[3];
rz(-2.9886768) q[3];
sx q[3];
rz(2.5482938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19330567) q[2];
sx q[2];
rz(-2.327658) q[2];
sx q[2];
rz(-1.8602271) q[2];
rz(-0.60503259) q[3];
sx q[3];
rz(-1.0993967) q[3];
sx q[3];
rz(1.9090778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.879377) q[0];
sx q[0];
rz(0.61304098) q[0];
rz(0.68060654) q[1];
sx q[1];
rz(-2.0304408) q[1];
sx q[1];
rz(-1.8162762) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3191736) q[0];
sx q[0];
rz(-2.1129916) q[0];
sx q[0];
rz(1.2577357) q[0];
rz(-pi) q[1];
rz(-1.4875796) q[2];
sx q[2];
rz(-1.7832547) q[2];
sx q[2];
rz(-2.0709289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5966353) q[1];
sx q[1];
rz(-1.4442634) q[1];
sx q[1];
rz(-0.56565779) q[1];
rz(-pi) q[2];
rz(0.57215055) q[3];
sx q[3];
rz(-0.66416262) q[3];
sx q[3];
rz(-0.45156878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1943872) q[2];
sx q[2];
rz(-1.2014061) q[2];
sx q[2];
rz(-2.3335333) q[2];
rz(0.64588532) q[3];
sx q[3];
rz(-0.26791993) q[3];
sx q[3];
rz(2.8382235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5658257) q[0];
sx q[0];
rz(-0.72661) q[0];
sx q[0];
rz(2.6988244) q[0];
rz(-0.59204656) q[1];
sx q[1];
rz(-2.4202085) q[1];
sx q[1];
rz(-0.17568849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4018049) q[0];
sx q[0];
rz(-1.8462088) q[0];
sx q[0];
rz(2.7729183) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4515516) q[2];
sx q[2];
rz(-1.6066666) q[2];
sx q[2];
rz(-0.6605688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86398839) q[1];
sx q[1];
rz(-0.93011198) q[1];
sx q[1];
rz(-1.5169657) q[1];
x q[2];
rz(-0.91276987) q[3];
sx q[3];
rz(-2.4173749) q[3];
sx q[3];
rz(-2.4778544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5438133) q[2];
sx q[2];
rz(-2.3976349) q[2];
sx q[2];
rz(2.234484) q[2];
rz(-0.90726566) q[3];
sx q[3];
rz(-1.3943358) q[3];
sx q[3];
rz(-3.0239014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1842781) q[0];
sx q[0];
rz(-0.90183455) q[0];
sx q[0];
rz(-0.19004518) q[0];
rz(1.5589145) q[1];
sx q[1];
rz(-1.5946486) q[1];
sx q[1];
rz(1.83439) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5763278) q[0];
sx q[0];
rz(-1.7143664) q[0];
sx q[0];
rz(1.2031892) q[0];
rz(-pi) q[1];
rz(2.2707269) q[2];
sx q[2];
rz(-1.3076837) q[2];
sx q[2];
rz(1.7990636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7293882) q[1];
sx q[1];
rz(-1.6784759) q[1];
sx q[1];
rz(3.1086224) q[1];
rz(-pi) q[2];
rz(-0.2449068) q[3];
sx q[3];
rz(-2.6953813) q[3];
sx q[3];
rz(3.0065766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.068453161) q[2];
sx q[2];
rz(-2.0544724) q[2];
sx q[2];
rz(2.9404409) q[2];
rz(-1.0226095) q[3];
sx q[3];
rz(-0.68251959) q[3];
sx q[3];
rz(1.539591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1744743) q[0];
sx q[0];
rz(-1.0459463) q[0];
sx q[0];
rz(-2.2651267) q[0];
rz(2.5190952) q[1];
sx q[1];
rz(-2.9298156) q[1];
sx q[1];
rz(2.7988787) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22607813) q[0];
sx q[0];
rz(-0.63185872) q[0];
sx q[0];
rz(1.4591239) q[0];
rz(-pi) q[1];
rz(1.9079952) q[2];
sx q[2];
rz(-2.7863481) q[2];
sx q[2];
rz(0.97081414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49149281) q[1];
sx q[1];
rz(-2.0180185) q[1];
sx q[1];
rz(-0.69112055) q[1];
rz(-pi) q[2];
rz(-0.56365074) q[3];
sx q[3];
rz(-0.4500126) q[3];
sx q[3];
rz(-1.8720418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9524625) q[2];
sx q[2];
rz(-1.6996982) q[2];
sx q[2];
rz(-2.6978037) q[2];
rz(-1.1187547) q[3];
sx q[3];
rz(-1.5951364) q[3];
sx q[3];
rz(0.55383033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8583869) q[0];
sx q[0];
rz(-1.5571742) q[0];
sx q[0];
rz(-1.5207186) q[0];
rz(-3.0781147) q[1];
sx q[1];
rz(-1.3451481) q[1];
sx q[1];
rz(1.4048911) q[1];
rz(-3.079698) q[2];
sx q[2];
rz(-1.9006922) q[2];
sx q[2];
rz(-3.047826) q[2];
rz(0.24334454) q[3];
sx q[3];
rz(-2.1235597) q[3];
sx q[3];
rz(-1.7490993) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
