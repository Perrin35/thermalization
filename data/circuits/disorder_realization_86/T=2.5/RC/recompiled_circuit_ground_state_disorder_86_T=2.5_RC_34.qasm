OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.13088432) q[0];
sx q[0];
rz(-0.67357981) q[0];
sx q[0];
rz(2.3186865) q[0];
rz(-1.6910488) q[1];
sx q[1];
rz(-2.4658642) q[1];
sx q[1];
rz(2.9339209) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7340644) q[0];
sx q[0];
rz(-0.82053608) q[0];
sx q[0];
rz(-1.090901) q[0];
x q[1];
rz(-0.34556371) q[2];
sx q[2];
rz(-1.0394382) q[2];
sx q[2];
rz(3.1329581) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.365578) q[1];
sx q[1];
rz(-1.509359) q[1];
sx q[1];
rz(0.12428026) q[1];
x q[2];
rz(-1.200233) q[3];
sx q[3];
rz(-1.5041122) q[3];
sx q[3];
rz(2.5197864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1755918) q[2];
sx q[2];
rz(-1.5755743) q[2];
sx q[2];
rz(1.9123745) q[2];
rz(-1.2980596) q[3];
sx q[3];
rz(-1.7652054) q[3];
sx q[3];
rz(-1.957533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9640279) q[0];
sx q[0];
rz(-0.25968817) q[0];
sx q[0];
rz(-2.2143256) q[0];
rz(-2.941653) q[1];
sx q[1];
rz(-1.6688469) q[1];
sx q[1];
rz(-0.98958579) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.171281) q[0];
sx q[0];
rz(-1.6785445) q[0];
sx q[0];
rz(1.2077483) q[0];
rz(-pi) q[1];
rz(-0.43439718) q[2];
sx q[2];
rz(-1.531562) q[2];
sx q[2];
rz(2.5920282) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11299388) q[1];
sx q[1];
rz(-1.6474287) q[1];
sx q[1];
rz(-3.0930685) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85736147) q[3];
sx q[3];
rz(-2.0318084) q[3];
sx q[3];
rz(-1.1198514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9280615) q[2];
sx q[2];
rz(-1.4175043) q[2];
sx q[2];
rz(0.35378635) q[2];
rz(2.1900322) q[3];
sx q[3];
rz(-0.74717251) q[3];
sx q[3];
rz(2.814754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3554409) q[0];
sx q[0];
rz(-0.94796258) q[0];
sx q[0];
rz(2.1308664) q[0];
rz(3.082869) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-0.40930632) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8008645) q[0];
sx q[0];
rz(-1.1363863) q[0];
sx q[0];
rz(-0.77462642) q[0];
x q[1];
rz(2.9761613) q[2];
sx q[2];
rz(-1.1647875) q[2];
sx q[2];
rz(2.5819786) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.18913705) q[1];
sx q[1];
rz(-2.5591095) q[1];
sx q[1];
rz(-2.7791695) q[1];
x q[2];
rz(2.6326551) q[3];
sx q[3];
rz(-1.5526315) q[3];
sx q[3];
rz(3.09336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.056444082) q[2];
sx q[2];
rz(-1.8935545) q[2];
sx q[2];
rz(0.15963456) q[2];
rz(-1.8917482) q[3];
sx q[3];
rz(-1.7227453) q[3];
sx q[3];
rz(1.0861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98642629) q[0];
sx q[0];
rz(-2.7041589) q[0];
sx q[0];
rz(-1.0945818) q[0];
rz(-1.9704341) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(-1.0135244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0243624) q[0];
sx q[0];
rz(-1.5268396) q[0];
sx q[0];
rz(-1.5323601) q[0];
rz(1.5640352) q[2];
sx q[2];
rz(-0.49443118) q[2];
sx q[2];
rz(1.0552366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78270528) q[1];
sx q[1];
rz(-0.26391477) q[1];
sx q[1];
rz(0.15586075) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0059399) q[3];
sx q[3];
rz(-1.6683104) q[3];
sx q[3];
rz(1.9103736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1363498) q[2];
sx q[2];
rz(-1.3815657) q[2];
sx q[2];
rz(1.3137777) q[2];
rz(-2.2037196) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(-3.042799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90792847) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(-1.044957) q[0];
rz(0.16432556) q[1];
sx q[1];
rz(-1.3427837) q[1];
sx q[1];
rz(-0.16990653) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92238802) q[0];
sx q[0];
rz(-1.6419852) q[0];
sx q[0];
rz(-1.3842596) q[0];
rz(-pi) q[1];
rz(-0.33268945) q[2];
sx q[2];
rz(-1.2342808) q[2];
sx q[2];
rz(2.3115186) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2881238) q[1];
sx q[1];
rz(-0.92902029) q[1];
sx q[1];
rz(0.30499129) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4154334) q[3];
sx q[3];
rz(-1.4724331) q[3];
sx q[3];
rz(-1.6189599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.56875151) q[2];
sx q[2];
rz(-0.70256394) q[2];
sx q[2];
rz(2.11002) q[2];
rz(-1.0860363) q[3];
sx q[3];
rz(-0.7998172) q[3];
sx q[3];
rz(-1.4097376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3373435) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(-1.7403437) q[0];
rz(0.42964545) q[1];
sx q[1];
rz(-0.88637561) q[1];
sx q[1];
rz(1.3053798) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8891762) q[0];
sx q[0];
rz(-2.2866797) q[0];
sx q[0];
rz(1.1578015) q[0];
rz(1.5611224) q[2];
sx q[2];
rz(-0.43894437) q[2];
sx q[2];
rz(-0.87359554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77338058) q[1];
sx q[1];
rz(-1.6016593) q[1];
sx q[1];
rz(0.24847757) q[1];
rz(1.8135728) q[3];
sx q[3];
rz(-1.7194028) q[3];
sx q[3];
rz(-0.44156238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0171011) q[2];
sx q[2];
rz(-1.93511) q[2];
sx q[2];
rz(2.1567832) q[2];
rz(-1.5752327) q[3];
sx q[3];
rz(-1.5519578) q[3];
sx q[3];
rz(-2.8563833) q[3];
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
rz(1.259909) q[0];
sx q[0];
rz(-2.8398828) q[0];
sx q[0];
rz(-1.786422) q[0];
rz(-3.1175218) q[1];
sx q[1];
rz(-1.6203974) q[1];
sx q[1];
rz(0.37758652) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51272362) q[0];
sx q[0];
rz(-2.2170904) q[0];
sx q[0];
rz(-1.3747526) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7837011) q[2];
sx q[2];
rz(-1.7434374) q[2];
sx q[2];
rz(1.5162692) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42492577) q[1];
sx q[1];
rz(-1.1831938) q[1];
sx q[1];
rz(2.7295154) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2144775) q[3];
sx q[3];
rz(-0.52442951) q[3];
sx q[3];
rz(2.7853109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54921237) q[2];
sx q[2];
rz(-0.33752957) q[2];
sx q[2];
rz(2.737992) q[2];
rz(0.18925439) q[3];
sx q[3];
rz(-1.0523825) q[3];
sx q[3];
rz(0.45690593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6117578) q[0];
sx q[0];
rz(-0.54946041) q[0];
sx q[0];
rz(0.31103617) q[0];
rz(-2.1850736) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(-0.32435736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3295591) q[0];
sx q[0];
rz(-2.5987465) q[0];
sx q[0];
rz(-0.37522786) q[0];
x q[1];
rz(-1.6300522) q[2];
sx q[2];
rz(-1.0612592) q[2];
sx q[2];
rz(-1.6409724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1970586) q[1];
sx q[1];
rz(-2.1558216) q[1];
sx q[1];
rz(2.6385175) q[1];
x q[2];
rz(0.37741669) q[3];
sx q[3];
rz(-2.2126865) q[3];
sx q[3];
rz(-2.9270594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4057464) q[2];
sx q[2];
rz(-2.3509071) q[2];
sx q[2];
rz(-0.22845593) q[2];
rz(-2.6089846) q[3];
sx q[3];
rz(-2.3753128) q[3];
sx q[3];
rz(2.2920091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.547895) q[0];
sx q[0];
rz(-2.6314647) q[0];
sx q[0];
rz(1.2702031) q[0];
rz(-2.5529329) q[1];
sx q[1];
rz(-1.207186) q[1];
sx q[1];
rz(-2.7587845) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4335121) q[0];
sx q[0];
rz(-1.5956586) q[0];
sx q[0];
rz(1.3415501) q[0];
x q[1];
rz(-1.8201939) q[2];
sx q[2];
rz(-1.1487085) q[2];
sx q[2];
rz(2.7673134) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11483604) q[1];
sx q[1];
rz(-0.72163218) q[1];
sx q[1];
rz(-2.9507157) q[1];
x q[2];
rz(0.270003) q[3];
sx q[3];
rz(-2.0473891) q[3];
sx q[3];
rz(1.6635513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1560912) q[2];
sx q[2];
rz(-1.5727377) q[2];
sx q[2];
rz(2.7413979) q[2];
rz(-0.24724809) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(2.7483773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.8409214) q[0];
sx q[0];
rz(-1.3341757) q[0];
sx q[0];
rz(-2.1260496) q[0];
rz(0.66889846) q[1];
sx q[1];
rz(-1.6872419) q[1];
sx q[1];
rz(1.5060172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43011623) q[0];
sx q[0];
rz(-1.073494) q[0];
sx q[0];
rz(2.3946752) q[0];
x q[1];
rz(-1.1349384) q[2];
sx q[2];
rz(-2.0761937) q[2];
sx q[2];
rz(-1.7541898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1387846) q[1];
sx q[1];
rz(-2.6600533) q[1];
sx q[1];
rz(-1.6549631) q[1];
rz(-pi) q[2];
rz(-0.90088441) q[3];
sx q[3];
rz(-0.79798079) q[3];
sx q[3];
rz(-0.98679286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16613913) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(-1.6157185) q[2];
rz(-2.8499991) q[3];
sx q[3];
rz(-0.44277954) q[3];
sx q[3];
rz(-2.7899138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780554) q[0];
sx q[0];
rz(-2.1778477) q[0];
sx q[0];
rz(0.5974593) q[0];
rz(0.28868227) q[1];
sx q[1];
rz(-0.79816993) q[1];
sx q[1];
rz(0.023963902) q[1];
rz(-1.2600675) q[2];
sx q[2];
rz(-1.1790397) q[2];
sx q[2];
rz(-2.9809748) q[2];
rz(2.6633429) q[3];
sx q[3];
rz(-0.50898715) q[3];
sx q[3];
rz(-0.57686808) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
