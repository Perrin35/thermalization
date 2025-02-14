OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49837056) q[0];
sx q[0];
rz(-1.4747488) q[0];
sx q[0];
rz(-0.21292444) q[0];
rz(-2.9990745) q[1];
sx q[1];
rz(-0.63894874) q[1];
sx q[1];
rz(0.45165935) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0395767) q[0];
sx q[0];
rz(-2.4964818) q[0];
sx q[0];
rz(-0.74104423) q[0];
rz(-pi) q[1];
rz(-0.84824382) q[2];
sx q[2];
rz(-2.179716) q[2];
sx q[2];
rz(-1.7681846) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38714007) q[1];
sx q[1];
rz(-1.3562227) q[1];
sx q[1];
rz(1.3438061) q[1];
x q[2];
rz(-2.0805243) q[3];
sx q[3];
rz(-0.76610111) q[3];
sx q[3];
rz(-1.776772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7529922) q[2];
sx q[2];
rz(-1.4326606) q[2];
sx q[2];
rz(0.60423744) q[2];
rz(-0.21449098) q[3];
sx q[3];
rz(-0.70591226) q[3];
sx q[3];
rz(0.96116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251505) q[0];
sx q[0];
rz(-2.3264628) q[0];
sx q[0];
rz(2.2863638) q[0];
rz(-1.1280355) q[1];
sx q[1];
rz(-0.74888343) q[1];
sx q[1];
rz(-1.523472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0013855502) q[0];
sx q[0];
rz(-1.1104776) q[0];
sx q[0];
rz(-2.1859474) q[0];
rz(-0.56407137) q[2];
sx q[2];
rz(-2.2937932) q[2];
sx q[2];
rz(-1.3644219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5828965) q[1];
sx q[1];
rz(-1.9882047) q[1];
sx q[1];
rz(-1.2701734) q[1];
rz(-pi) q[2];
rz(0.43782708) q[3];
sx q[3];
rz(-1.722933) q[3];
sx q[3];
rz(0.9532387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7889859) q[2];
sx q[2];
rz(-1.8800198) q[2];
sx q[2];
rz(-1.2471586) q[2];
rz(-2.9712307) q[3];
sx q[3];
rz(-2.7046461) q[3];
sx q[3];
rz(-0.40444571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929297) q[0];
sx q[0];
rz(-0.12032838) q[0];
sx q[0];
rz(2.2509101) q[0];
rz(1.0825253) q[1];
sx q[1];
rz(-2.4286916) q[1];
sx q[1];
rz(-3.069186) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88806498) q[0];
sx q[0];
rz(-2.088659) q[0];
sx q[0];
rz(-2.6648471) q[0];
rz(-pi) q[1];
rz(2.7909325) q[2];
sx q[2];
rz(-1.114801) q[2];
sx q[2];
rz(0.58619546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9074633) q[1];
sx q[1];
rz(-1.0765706) q[1];
sx q[1];
rz(2.0496561) q[1];
rz(1.6909353) q[3];
sx q[3];
rz(-1.9124219) q[3];
sx q[3];
rz(1.6610314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.801164) q[2];
sx q[2];
rz(-1.9457996) q[2];
sx q[2];
rz(2.8677531) q[2];
rz(1.0770816) q[3];
sx q[3];
rz(-1.8718953) q[3];
sx q[3];
rz(-2.4382408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3649789) q[0];
sx q[0];
rz(-0.7319428) q[0];
sx q[0];
rz(-1.8044385) q[0];
rz(2.8058167) q[1];
sx q[1];
rz(-1.8700721) q[1];
sx q[1];
rz(1.5912067) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5364658) q[0];
sx q[0];
rz(-0.7881645) q[0];
sx q[0];
rz(0.76351662) q[0];
x q[1];
rz(-2.6350722) q[2];
sx q[2];
rz(-1.6088076) q[2];
sx q[2];
rz(2.4354629) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.76564255) q[1];
sx q[1];
rz(-0.11184622) q[1];
sx q[1];
rz(1.8814398) q[1];
rz(2.4404686) q[3];
sx q[3];
rz(-1.935734) q[3];
sx q[3];
rz(-2.622245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7871899) q[2];
sx q[2];
rz(-1.411011) q[2];
sx q[2];
rz(-1.2350941) q[2];
rz(0.45803329) q[3];
sx q[3];
rz(-2.1221275) q[3];
sx q[3];
rz(0.73634806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6422727) q[0];
sx q[0];
rz(-0.90773931) q[0];
sx q[0];
rz(0.15550144) q[0];
rz(2.3606965) q[1];
sx q[1];
rz(-0.28762329) q[1];
sx q[1];
rz(0.23316613) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6885666) q[0];
sx q[0];
rz(-0.84491623) q[0];
sx q[0];
rz(0.19832439) q[0];
x q[1];
rz(0.16712971) q[2];
sx q[2];
rz(-2.4592113) q[2];
sx q[2];
rz(-0.025957195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.609954) q[1];
sx q[1];
rz(-1.8401658) q[1];
sx q[1];
rz(-0.1166719) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.030117) q[3];
sx q[3];
rz(-1.0305163) q[3];
sx q[3];
rz(-0.92898166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1843725) q[2];
sx q[2];
rz(-1.5868712) q[2];
sx q[2];
rz(0.30259821) q[2];
rz(0.042081632) q[3];
sx q[3];
rz(-2.849597) q[3];
sx q[3];
rz(-0.80872768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0724723) q[0];
sx q[0];
rz(-1.1061677) q[0];
sx q[0];
rz(0.99232596) q[0];
rz(0.26556695) q[1];
sx q[1];
rz(-2.3873603) q[1];
sx q[1];
rz(0.011836424) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3385394) q[0];
sx q[0];
rz(-1.526471) q[0];
sx q[0];
rz(-1.2903777) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1407848) q[2];
sx q[2];
rz(-1.3234183) q[2];
sx q[2];
rz(0.22817366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1975147) q[1];
sx q[1];
rz(-1.6374365) q[1];
sx q[1];
rz(1.2032894) q[1];
rz(-pi) q[2];
rz(1.2740259) q[3];
sx q[3];
rz(-2.9057876) q[3];
sx q[3];
rz(0.91487003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5609694) q[2];
sx q[2];
rz(-0.89726609) q[2];
sx q[2];
rz(0.90274367) q[2];
rz(0.48815253) q[3];
sx q[3];
rz(-1.3267696) q[3];
sx q[3];
rz(1.7972402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1667267) q[0];
sx q[0];
rz(-2.6868197) q[0];
sx q[0];
rz(1.3653261) q[0];
rz(2.4873554) q[1];
sx q[1];
rz(-2.3932494) q[1];
sx q[1];
rz(-1.0985589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7099534) q[0];
sx q[0];
rz(-1.1134143) q[0];
sx q[0];
rz(-2.4541992) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8945872) q[2];
sx q[2];
rz(-1.5135458) q[2];
sx q[2];
rz(1.0047785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85118402) q[1];
sx q[1];
rz(-0.99857226) q[1];
sx q[1];
rz(1.6838724) q[1];
rz(-pi) q[2];
rz(2.8756254) q[3];
sx q[3];
rz(-2.1823357) q[3];
sx q[3];
rz(-2.3872294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4174623) q[2];
sx q[2];
rz(-1.7899568) q[2];
sx q[2];
rz(-2.2473118) q[2];
rz(1.2096679) q[3];
sx q[3];
rz(-3.0318048) q[3];
sx q[3];
rz(-2.3144058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9342448) q[0];
sx q[0];
rz(-1.8266015) q[0];
sx q[0];
rz(-0.19499245) q[0];
rz(0.64942819) q[1];
sx q[1];
rz(-0.94562999) q[1];
sx q[1];
rz(1.6920998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6882807) q[0];
sx q[0];
rz(-1.9460228) q[0];
sx q[0];
rz(-0.58102258) q[0];
rz(-pi) q[1];
rz(-0.56153293) q[2];
sx q[2];
rz(-0.61164633) q[2];
sx q[2];
rz(1.2831519) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8713432) q[1];
sx q[1];
rz(-0.91971524) q[1];
sx q[1];
rz(3.0581942) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4420801) q[3];
sx q[3];
rz(-1.2815003) q[3];
sx q[3];
rz(-2.4955179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.080716915) q[2];
sx q[2];
rz(-0.94451153) q[2];
sx q[2];
rz(-0.057849217) q[2];
rz(3.0625693) q[3];
sx q[3];
rz(-1.9130324) q[3];
sx q[3];
rz(1.2270989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77646112) q[0];
sx q[0];
rz(-1.6420028) q[0];
sx q[0];
rz(2.814433) q[0];
rz(-2.6311334) q[1];
sx q[1];
rz(-1.8233428) q[1];
sx q[1];
rz(-3.0697451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15639977) q[0];
sx q[0];
rz(-1.4738393) q[0];
sx q[0];
rz(-0.18228874) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4445093) q[2];
sx q[2];
rz(-1.0490745) q[2];
sx q[2];
rz(1.9283805) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0409826) q[1];
sx q[1];
rz(-2.1400053) q[1];
sx q[1];
rz(2.3952574) q[1];
rz(-0.44485299) q[3];
sx q[3];
rz(-2.3844845) q[3];
sx q[3];
rz(-3.0187154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0923126) q[2];
sx q[2];
rz(-1.6593554) q[2];
sx q[2];
rz(1.9149038) q[2];
rz(-0.53931326) q[3];
sx q[3];
rz(-1.156811) q[3];
sx q[3];
rz(1.4815319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6698089) q[0];
sx q[0];
rz(-2.0150549) q[0];
sx q[0];
rz(-0.54927611) q[0];
rz(-1.1726941) q[1];
sx q[1];
rz(-1.6212308) q[1];
sx q[1];
rz(-1.2991914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8912997) q[0];
sx q[0];
rz(-1.5101103) q[0];
sx q[0];
rz(2.9461198) q[0];
x q[1];
rz(0.46165919) q[2];
sx q[2];
rz(-0.73349059) q[2];
sx q[2];
rz(-1.5635827) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.047094237) q[1];
sx q[1];
rz(-2.9278122) q[1];
sx q[1];
rz(-0.70014145) q[1];
rz(-pi) q[2];
rz(2.5985093) q[3];
sx q[3];
rz(-1.833931) q[3];
sx q[3];
rz(-3.0656726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4486763) q[2];
sx q[2];
rz(-2.0557949) q[2];
sx q[2];
rz(2.9936301) q[2];
rz(-0.97370094) q[3];
sx q[3];
rz(-0.86108834) q[3];
sx q[3];
rz(-2.2073943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14514087) q[0];
sx q[0];
rz(-1.8358163) q[0];
sx q[0];
rz(-1.3124574) q[0];
rz(-1.8295857) q[1];
sx q[1];
rz(-0.94999718) q[1];
sx q[1];
rz(-2.2012262) q[1];
rz(-0.13036556) q[2];
sx q[2];
rz(-1.833537) q[2];
sx q[2];
rz(-0.037787211) q[2];
rz(0.14518006) q[3];
sx q[3];
rz(-0.96031453) q[3];
sx q[3];
rz(1.7688383) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
