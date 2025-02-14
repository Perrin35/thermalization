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
rz(2.5026439) q[1];
sx q[1];
rz(12.114711) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0395767) q[0];
sx q[0];
rz(-0.64511089) q[0];
sx q[0];
rz(2.4005484) q[0];
rz(-pi) q[1];
rz(2.2933488) q[2];
sx q[2];
rz(-0.96187667) q[2];
sx q[2];
rz(-1.3734081) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2131552) q[1];
sx q[1];
rz(-0.31107956) q[1];
sx q[1];
rz(-2.3401287) q[1];
rz(-pi) q[2];
rz(0.87224023) q[3];
sx q[3];
rz(-1.2256825) q[3];
sx q[3];
rz(0.17696571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3886005) q[2];
sx q[2];
rz(-1.4326606) q[2];
sx q[2];
rz(-2.5373552) q[2];
rz(-0.21449098) q[3];
sx q[3];
rz(-2.4356804) q[3];
sx q[3];
rz(-0.96116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.71644217) q[0];
sx q[0];
rz(-0.81512988) q[0];
sx q[0];
rz(-2.2863638) q[0];
rz(-2.0135571) q[1];
sx q[1];
rz(-2.3927092) q[1];
sx q[1];
rz(-1.523472) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0013855502) q[0];
sx q[0];
rz(-2.0311151) q[0];
sx q[0];
rz(-2.1859474) q[0];
rz(-pi) q[1];
rz(-0.56407137) q[2];
sx q[2];
rz(-2.2937932) q[2];
sx q[2];
rz(-1.3644219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.93003004) q[1];
sx q[1];
rz(-0.50917823) q[1];
sx q[1];
rz(-2.5528564) q[1];
x q[2];
rz(-0.43782708) q[3];
sx q[3];
rz(-1.722933) q[3];
sx q[3];
rz(2.188354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7889859) q[2];
sx q[2];
rz(-1.2615729) q[2];
sx q[2];
rz(1.2471586) q[2];
rz(-0.17036197) q[3];
sx q[3];
rz(-0.43694654) q[3];
sx q[3];
rz(-0.40444571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4929297) q[0];
sx q[0];
rz(-0.12032838) q[0];
sx q[0];
rz(0.89068252) q[0];
rz(1.0825253) q[1];
sx q[1];
rz(-0.71290103) q[1];
sx q[1];
rz(-0.072406702) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601357) q[0];
sx q[0];
rz(-2.4527672) q[0];
sx q[0];
rz(-2.248855) q[0];
rz(-pi) q[1];
rz(2.7909325) q[2];
sx q[2];
rz(-1.114801) q[2];
sx q[2];
rz(-2.5553972) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9074633) q[1];
sx q[1];
rz(-1.0765706) q[1];
sx q[1];
rz(-2.0496561) q[1];
rz(-pi) q[2];
rz(-0.32510783) q[3];
sx q[3];
rz(-2.7802411) q[3];
sx q[3];
rz(2.0068866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.801164) q[2];
sx q[2];
rz(-1.195793) q[2];
sx q[2];
rz(0.27383956) q[2];
rz(1.0770816) q[3];
sx q[3];
rz(-1.8718953) q[3];
sx q[3];
rz(0.70335189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7766137) q[0];
sx q[0];
rz(-0.7319428) q[0];
sx q[0];
rz(1.3371542) q[0];
rz(-2.8058167) q[1];
sx q[1];
rz(-1.2715205) q[1];
sx q[1];
rz(1.5912067) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44061892) q[0];
sx q[0];
rz(-1.0583726) q[0];
sx q[0];
rz(-2.5133564) q[0];
rz(-2.6350722) q[2];
sx q[2];
rz(-1.6088076) q[2];
sx q[2];
rz(2.4354629) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.027615) q[1];
sx q[1];
rz(-1.5366728) q[1];
sx q[1];
rz(1.464262) q[1];
x q[2];
rz(0.70112409) q[3];
sx q[3];
rz(-1.935734) q[3];
sx q[3];
rz(-0.51934767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35440272) q[2];
sx q[2];
rz(-1.7305817) q[2];
sx q[2];
rz(-1.9064986) q[2];
rz(-0.45803329) q[3];
sx q[3];
rz(-1.0194651) q[3];
sx q[3];
rz(-2.4052446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4993199) q[0];
sx q[0];
rz(-2.2338533) q[0];
sx q[0];
rz(-0.15550144) q[0];
rz(2.3606965) q[1];
sx q[1];
rz(-2.8539694) q[1];
sx q[1];
rz(2.9084265) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6885666) q[0];
sx q[0];
rz(-0.84491623) q[0];
sx q[0];
rz(-2.9432683) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16712971) q[2];
sx q[2];
rz(-0.68238133) q[2];
sx q[2];
rz(-3.1156355) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0079769) q[1];
sx q[1];
rz(-1.6832428) q[1];
sx q[1];
rz(-1.2996718) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70919819) q[3];
sx q[3];
rz(-0.74477421) q[3];
sx q[3];
rz(1.7913558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1843725) q[2];
sx q[2];
rz(-1.5547215) q[2];
sx q[2];
rz(2.8389944) q[2];
rz(-0.042081632) q[3];
sx q[3];
rz(-2.849597) q[3];
sx q[3];
rz(0.80872768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0724723) q[0];
sx q[0];
rz(-1.1061677) q[0];
sx q[0];
rz(0.99232596) q[0];
rz(-0.26556695) q[1];
sx q[1];
rz(-0.75423232) q[1];
sx q[1];
rz(0.011836424) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3385394) q[0];
sx q[0];
rz(-1.526471) q[0];
sx q[0];
rz(1.851215) q[0];
rz(-1.0008079) q[2];
sx q[2];
rz(-1.8181744) q[2];
sx q[2];
rz(2.913419) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34764987) q[1];
sx q[1];
rz(-1.2041438) q[1];
sx q[1];
rz(-0.071392752) q[1];
x q[2];
rz(-1.2740259) q[3];
sx q[3];
rz(-0.23580509) q[3];
sx q[3];
rz(0.91487003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5806233) q[2];
sx q[2];
rz(-2.2443266) q[2];
sx q[2];
rz(-0.90274367) q[2];
rz(-2.6534401) q[3];
sx q[3];
rz(-1.3267696) q[3];
sx q[3];
rz(1.7972402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.97486597) q[0];
sx q[0];
rz(-0.45477295) q[0];
sx q[0];
rz(-1.7762666) q[0];
rz(2.4873554) q[1];
sx q[1];
rz(-0.74834329) q[1];
sx q[1];
rz(1.0985589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43163928) q[0];
sx q[0];
rz(-1.1134143) q[0];
sx q[0];
rz(2.4541992) q[0];
rz(-1.5117579) q[2];
sx q[2];
rz(-1.8173886) q[2];
sx q[2];
rz(0.55159071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0579018) q[1];
sx q[1];
rz(-2.5595287) q[1];
sx q[1];
rz(2.9681724) q[1];
x q[2];
rz(0.94233108) q[3];
sx q[3];
rz(-1.7876995) q[3];
sx q[3];
rz(-2.1700117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4174623) q[2];
sx q[2];
rz(-1.7899568) q[2];
sx q[2];
rz(-0.89428085) q[2];
rz(1.9319247) q[3];
sx q[3];
rz(-3.0318048) q[3];
sx q[3];
rz(2.3144058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2073479) q[0];
sx q[0];
rz(-1.3149911) q[0];
sx q[0];
rz(-0.19499245) q[0];
rz(0.64942819) q[1];
sx q[1];
rz(-0.94562999) q[1];
sx q[1];
rz(1.6920998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7879542) q[0];
sx q[0];
rz(-1.0348085) q[0];
sx q[0];
rz(1.130442) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5800597) q[2];
sx q[2];
rz(-2.5299463) q[2];
sx q[2];
rz(1.2831519) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2499293) q[1];
sx q[1];
rz(-1.504487) q[1];
sx q[1];
rz(-2.2235566) q[1];
rz(-0.6995126) q[3];
sx q[3];
rz(-1.2815003) q[3];
sx q[3];
rz(-2.4955179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.080716915) q[2];
sx q[2];
rz(-0.94451153) q[2];
sx q[2];
rz(-3.0837434) q[2];
rz(3.0625693) q[3];
sx q[3];
rz(-1.2285602) q[3];
sx q[3];
rz(1.9144937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3651315) q[0];
sx q[0];
rz(-1.6420028) q[0];
sx q[0];
rz(-0.32715964) q[0];
rz(2.6311334) q[1];
sx q[1];
rz(-1.8233428) q[1];
sx q[1];
rz(3.0697451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8979683) q[0];
sx q[0];
rz(-0.20621696) q[0];
sx q[0];
rz(-2.6491524) q[0];
rz(-pi) q[1];
rz(-2.2141404) q[2];
sx q[2];
rz(-0.98053741) q[2];
sx q[2];
rz(3.1038499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0409826) q[1];
sx q[1];
rz(-2.1400053) q[1];
sx q[1];
rz(-0.74633523) q[1];
x q[2];
rz(1.9570146) q[3];
sx q[3];
rz(-0.90208331) q[3];
sx q[3];
rz(0.45763256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0492801) q[2];
sx q[2];
rz(-1.6593554) q[2];
sx q[2];
rz(1.9149038) q[2];
rz(-2.6022794) q[3];
sx q[3];
rz(-1.9847816) q[3];
sx q[3];
rz(-1.6600608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4717838) q[0];
sx q[0];
rz(-2.0150549) q[0];
sx q[0];
rz(-0.54927611) q[0];
rz(1.9688985) q[1];
sx q[1];
rz(-1.5203618) q[1];
sx q[1];
rz(1.2991914) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25029297) q[0];
sx q[0];
rz(-1.5101103) q[0];
sx q[0];
rz(-2.9461198) q[0];
x q[1];
rz(-1.9525398) q[2];
sx q[2];
rz(-2.2135268) q[2];
sx q[2];
rz(2.1537154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.047094237) q[1];
sx q[1];
rz(-2.9278122) q[1];
sx q[1];
rz(2.4414512) q[1];
x q[2];
rz(0.54308335) q[3];
sx q[3];
rz(-1.3076617) q[3];
sx q[3];
rz(0.07592003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6929163) q[2];
sx q[2];
rz(-1.0857978) q[2];
sx q[2];
rz(-2.9936301) q[2];
rz(-0.97370094) q[3];
sx q[3];
rz(-0.86108834) q[3];
sx q[3];
rz(0.93419832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14514087) q[0];
sx q[0];
rz(-1.8358163) q[0];
sx q[0];
rz(-1.3124574) q[0];
rz(1.3120069) q[1];
sx q[1];
rz(-0.94999718) q[1];
sx q[1];
rz(-2.2012262) q[1];
rz(-0.13036556) q[2];
sx q[2];
rz(-1.833537) q[2];
sx q[2];
rz(-0.037787211) q[2];
rz(2.9964126) q[3];
sx q[3];
rz(-2.1812781) q[3];
sx q[3];
rz(-1.3727544) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
