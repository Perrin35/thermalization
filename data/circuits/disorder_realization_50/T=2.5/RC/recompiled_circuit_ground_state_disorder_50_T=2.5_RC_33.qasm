OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6432221) q[0];
sx q[0];
rz(-1.6668439) q[0];
sx q[0];
rz(-2.9286682) q[0];
rz(-2.9990745) q[1];
sx q[1];
rz(-0.63894874) q[1];
sx q[1];
rz(-2.6899333) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925326) q[0];
sx q[0];
rz(-1.1111759) q[0];
sx q[0];
rz(-1.100774) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3826936) q[2];
sx q[2];
rz(-2.2335608) q[2];
sx q[2];
rz(-2.3686705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7544526) q[1];
sx q[1];
rz(-1.78537) q[1];
sx q[1];
rz(-1.3438061) q[1];
rz(-pi) q[2];
rz(1.0610683) q[3];
sx q[3];
rz(-0.76610111) q[3];
sx q[3];
rz(1.3648206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7529922) q[2];
sx q[2];
rz(-1.4326606) q[2];
sx q[2];
rz(0.60423744) q[2];
rz(2.9271017) q[3];
sx q[3];
rz(-0.70591226) q[3];
sx q[3];
rz(0.96116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71644217) q[0];
sx q[0];
rz(-0.81512988) q[0];
sx q[0];
rz(-0.85522884) q[0];
rz(-2.0135571) q[1];
sx q[1];
rz(-0.74888343) q[1];
sx q[1];
rz(1.523472) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0013855502) q[0];
sx q[0];
rz(-2.0311151) q[0];
sx q[0];
rz(-2.1859474) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1155384) q[2];
sx q[2];
rz(-2.2570809) q[2];
sx q[2];
rz(1.0141745) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93003004) q[1];
sx q[1];
rz(-0.50917823) q[1];
sx q[1];
rz(0.58873621) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4030966) q[3];
sx q[3];
rz(-2.0032231) q[3];
sx q[3];
rz(-2.5948615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35260674) q[2];
sx q[2];
rz(-1.2615729) q[2];
sx q[2];
rz(1.2471586) q[2];
rz(2.9712307) q[3];
sx q[3];
rz(-2.7046461) q[3];
sx q[3];
rz(0.40444571) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64866292) q[0];
sx q[0];
rz(-3.0212643) q[0];
sx q[0];
rz(2.2509101) q[0];
rz(-1.0825253) q[1];
sx q[1];
rz(-2.4286916) q[1];
sx q[1];
rz(3.069186) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0601357) q[0];
sx q[0];
rz(-2.4527672) q[0];
sx q[0];
rz(-0.89273767) q[0];
rz(-pi) q[1];
rz(-0.35066013) q[2];
sx q[2];
rz(-2.0267916) q[2];
sx q[2];
rz(-0.58619546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.23412936) q[1];
sx q[1];
rz(-2.0650221) q[1];
sx q[1];
rz(1.0919366) q[1];
x q[2];
rz(-2.8164848) q[3];
sx q[3];
rz(-2.7802411) q[3];
sx q[3];
rz(-2.0068866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34042865) q[2];
sx q[2];
rz(-1.9457996) q[2];
sx q[2];
rz(0.27383956) q[2];
rz(-1.0770816) q[3];
sx q[3];
rz(-1.2696973) q[3];
sx q[3];
rz(-2.4382408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3649789) q[0];
sx q[0];
rz(-2.4096498) q[0];
sx q[0];
rz(1.3371542) q[0];
rz(-2.8058167) q[1];
sx q[1];
rz(-1.2715205) q[1];
sx q[1];
rz(1.5912067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6051269) q[0];
sx q[0];
rz(-2.3534282) q[0];
sx q[0];
rz(-0.76351662) q[0];
rz(-pi) q[1];
rz(1.6142593) q[2];
sx q[2];
rz(-1.0646766) q[2];
sx q[2];
rz(-2.2558457) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3759501) q[1];
sx q[1];
rz(-0.11184622) q[1];
sx q[1];
rz(1.2601529) q[1];
rz(2.4404686) q[3];
sx q[3];
rz(-1.935734) q[3];
sx q[3];
rz(0.51934767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7871899) q[2];
sx q[2];
rz(-1.7305817) q[2];
sx q[2];
rz(-1.9064986) q[2];
rz(-2.6835594) q[3];
sx q[3];
rz(-1.0194651) q[3];
sx q[3];
rz(2.4052446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6422727) q[0];
sx q[0];
rz(-0.90773931) q[0];
sx q[0];
rz(-0.15550144) q[0];
rz(2.3606965) q[1];
sx q[1];
rz(-0.28762329) q[1];
sx q[1];
rz(0.23316613) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45302603) q[0];
sx q[0];
rz(-2.2966764) q[0];
sx q[0];
rz(-2.9432683) q[0];
x q[1];
rz(-0.67552394) q[2];
sx q[2];
rz(-1.6758989) q[2];
sx q[2];
rz(-1.6750248) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1951016) q[1];
sx q[1];
rz(-0.29298654) q[1];
sx q[1];
rz(1.9698214) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.030117) q[3];
sx q[3];
rz(-2.1110764) q[3];
sx q[3];
rz(-2.212611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1843725) q[2];
sx q[2];
rz(-1.5868712) q[2];
sx q[2];
rz(-2.8389944) q[2];
rz(-0.042081632) q[3];
sx q[3];
rz(-2.849597) q[3];
sx q[3];
rz(-2.332865) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06912032) q[0];
sx q[0];
rz(-1.1061677) q[0];
sx q[0];
rz(2.1492667) q[0];
rz(2.8760257) q[1];
sx q[1];
rz(-2.3873603) q[1];
sx q[1];
rz(3.1297562) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079601057) q[0];
sx q[0];
rz(-2.8577836) q[0];
sx q[0];
rz(1.4118836) q[0];
x q[1];
rz(0.29143362) q[2];
sx q[2];
rz(-1.0202006) q[2];
sx q[2];
rz(1.498298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1975147) q[1];
sx q[1];
rz(-1.5041562) q[1];
sx q[1];
rz(1.9383033) q[1];
rz(-pi) q[2];
rz(-1.2740259) q[3];
sx q[3];
rz(-0.23580509) q[3];
sx q[3];
rz(0.91487003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5806233) q[2];
sx q[2];
rz(-0.89726609) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97486597) q[0];
sx q[0];
rz(-2.6868197) q[0];
sx q[0];
rz(-1.7762666) q[0];
rz(-2.4873554) q[1];
sx q[1];
rz(-0.74834329) q[1];
sx q[1];
rz(-1.0985589) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64564686) q[0];
sx q[0];
rz(-2.3370727) q[0];
sx q[0];
rz(-0.65973892) q[0];
rz(-pi) q[1];
rz(-1.6298348) q[2];
sx q[2];
rz(-1.8173886) q[2];
sx q[2];
rz(-0.55159071) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.65819627) q[1];
sx q[1];
rz(-1.475793) q[1];
sx q[1];
rz(0.57514455) q[1];
x q[2];
rz(-1.9294268) q[3];
sx q[3];
rz(-2.4815668) q[3];
sx q[3];
rz(2.8302659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72413033) q[2];
sx q[2];
rz(-1.7899568) q[2];
sx q[2];
rz(-2.2473118) q[2];
rz(-1.2096679) q[3];
sx q[3];
rz(-3.0318048) q[3];
sx q[3];
rz(-0.82718682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.2073479) q[0];
sx q[0];
rz(-1.3149911) q[0];
sx q[0];
rz(0.19499245) q[0];
rz(0.64942819) q[1];
sx q[1];
rz(-2.1959627) q[1];
sx q[1];
rz(1.4494928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45331199) q[0];
sx q[0];
rz(-1.1955698) q[0];
sx q[0];
rz(0.58102258) q[0];
rz(0.53575164) q[2];
sx q[2];
rz(-1.2600586) q[2];
sx q[2];
rz(-0.18795943) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7342775) q[1];
sx q[1];
rz(-2.4859634) q[1];
sx q[1];
rz(-1.6796965) q[1];
rz(0.43303135) q[3];
sx q[3];
rz(-2.3940475) q[3];
sx q[3];
rz(1.8900035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.080716915) q[2];
sx q[2];
rz(-2.1970811) q[2];
sx q[2];
rz(0.057849217) q[2];
rz(3.0625693) q[3];
sx q[3];
rz(-1.9130324) q[3];
sx q[3];
rz(1.2270989) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77646112) q[0];
sx q[0];
rz(-1.4995898) q[0];
sx q[0];
rz(-0.32715964) q[0];
rz(-2.6311334) q[1];
sx q[1];
rz(-1.8233428) q[1];
sx q[1];
rz(-3.0697451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15639977) q[0];
sx q[0];
rz(-1.4738393) q[0];
sx q[0];
rz(-0.18228874) q[0];
rz(-pi) q[1];
rz(-2.4113095) q[2];
sx q[2];
rz(-0.84368333) q[2];
sx q[2];
rz(-2.2471225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0409826) q[1];
sx q[1];
rz(-1.0015873) q[1];
sx q[1];
rz(-0.74633523) q[1];
x q[2];
rz(0.44485299) q[3];
sx q[3];
rz(-0.75710812) q[3];
sx q[3];
rz(0.12287724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0492801) q[2];
sx q[2];
rz(-1.6593554) q[2];
sx q[2];
rz(-1.9149038) q[2];
rz(2.6022794) q[3];
sx q[3];
rz(-1.9847816) q[3];
sx q[3];
rz(-1.4815319) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6698089) q[0];
sx q[0];
rz(-1.1265378) q[0];
sx q[0];
rz(2.5923165) q[0];
rz(1.1726941) q[1];
sx q[1];
rz(-1.6212308) q[1];
sx q[1];
rz(1.2991914) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25029297) q[0];
sx q[0];
rz(-1.5101103) q[0];
sx q[0];
rz(0.19547284) q[0];
rz(-1.1890529) q[2];
sx q[2];
rz(-0.92806584) q[2];
sx q[2];
rz(2.1537154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92907266) q[1];
sx q[1];
rz(-1.4336695) q[1];
sx q[1];
rz(-0.16452505) q[1];
rz(-pi) q[2];
rz(0.48051832) q[3];
sx q[3];
rz(-2.5439265) q[3];
sx q[3];
rz(1.0880119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6929163) q[2];
sx q[2];
rz(-2.0557949) q[2];
sx q[2];
rz(2.9936301) q[2];
rz(-2.1678917) q[3];
sx q[3];
rz(-0.86108834) q[3];
sx q[3];
rz(2.2073943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9964518) q[0];
sx q[0];
rz(-1.3057764) q[0];
sx q[0];
rz(1.8291352) q[0];
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
