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
rz(4.8084365) q[0];
sx q[0];
rz(9.2118535) q[0];
rz(-2.9990745) q[1];
sx q[1];
rz(-0.63894874) q[1];
sx q[1];
rz(0.45165935) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.102016) q[0];
sx q[0];
rz(-2.4964818) q[0];
sx q[0];
rz(2.4005484) q[0];
rz(-pi) q[1];
rz(0.75889905) q[2];
sx q[2];
rz(-0.90803185) q[2];
sx q[2];
rz(-0.77292216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38714007) q[1];
sx q[1];
rz(-1.78537) q[1];
sx q[1];
rz(1.7977865) q[1];
rz(2.0805243) q[3];
sx q[3];
rz(-0.76610111) q[3];
sx q[3];
rz(-1.3648206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.3886005) q[2];
sx q[2];
rz(-1.4326606) q[2];
sx q[2];
rz(-2.5373552) q[2];
rz(0.21449098) q[3];
sx q[3];
rz(-0.70591226) q[3];
sx q[3];
rz(2.1804325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71644217) q[0];
sx q[0];
rz(-2.3264628) q[0];
sx q[0];
rz(0.85522884) q[0];
rz(-2.0135571) q[1];
sx q[1];
rz(-0.74888343) q[1];
sx q[1];
rz(-1.6181207) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.010941) q[0];
sx q[0];
rz(-2.391577) q[0];
sx q[0];
rz(2.2806243) q[0];
rz(-2.3777826) q[2];
sx q[2];
rz(-1.1583185) q[2];
sx q[2];
rz(-0.19004059) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11292371) q[1];
sx q[1];
rz(-1.296685) q[1];
sx q[1];
rz(2.7069189) q[1];
rz(2.7945943) q[3];
sx q[3];
rz(-2.6797077) q[3];
sx q[3];
rz(-2.210946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35260674) q[2];
sx q[2];
rz(-1.2615729) q[2];
sx q[2];
rz(-1.2471586) q[2];
rz(2.9712307) q[3];
sx q[3];
rz(-2.7046461) q[3];
sx q[3];
rz(-2.7371469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929297) q[0];
sx q[0];
rz(-0.12032838) q[0];
sx q[0];
rz(2.2509101) q[0];
rz(2.0590674) q[1];
sx q[1];
rz(-0.71290103) q[1];
sx q[1];
rz(-3.069186) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081456979) q[0];
sx q[0];
rz(-0.68882548) q[0];
sx q[0];
rz(2.248855) q[0];
rz(-0.95981284) q[2];
sx q[2];
rz(-0.56761203) q[2];
sx q[2];
rz(-1.8622514) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5634656) q[1];
sx q[1];
rz(-1.153128) q[1];
sx q[1];
rz(-0.54563569) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-0.34042865) q[2];
sx q[2];
rz(-1.9457996) q[2];
sx q[2];
rz(2.8677531) q[2];
rz(-2.0645111) q[3];
sx q[3];
rz(-1.8718953) q[3];
sx q[3];
rz(-2.4382408) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7766137) q[0];
sx q[0];
rz(-0.7319428) q[0];
sx q[0];
rz(1.8044385) q[0];
rz(-0.33577597) q[1];
sx q[1];
rz(-1.2715205) q[1];
sx q[1];
rz(1.550386) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7009737) q[0];
sx q[0];
rz(-2.0832201) q[0];
sx q[0];
rz(-0.62823624) q[0];
rz(2.6350722) q[2];
sx q[2];
rz(-1.532785) q[2];
sx q[2];
rz(2.4354629) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76564255) q[1];
sx q[1];
rz(-0.11184622) q[1];
sx q[1];
rz(1.8814398) q[1];
rz(2.0344355) q[3];
sx q[3];
rz(-2.2175466) q[3];
sx q[3];
rz(-1.3440901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35440272) q[2];
sx q[2];
rz(-1.411011) q[2];
sx q[2];
rz(-1.2350941) q[2];
rz(-2.6835594) q[3];
sx q[3];
rz(-1.0194651) q[3];
sx q[3];
rz(-0.73634806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6422727) q[0];
sx q[0];
rz(-2.2338533) q[0];
sx q[0];
rz(0.15550144) q[0];
rz(-2.3606965) q[1];
sx q[1];
rz(-2.8539694) q[1];
sx q[1];
rz(-2.9084265) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45302603) q[0];
sx q[0];
rz(-0.84491623) q[0];
sx q[0];
rz(-2.9432683) q[0];
rz(-pi) q[1];
rz(2.4660687) q[2];
sx q[2];
rz(-1.4656938) q[2];
sx q[2];
rz(-1.4665678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53163869) q[1];
sx q[1];
rz(-1.8401658) q[1];
sx q[1];
rz(0.1166719) q[1];
rz(-2.5311337) q[3];
sx q[3];
rz(-1.1136331) q[3];
sx q[3];
rz(2.7992953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1843725) q[2];
sx q[2];
rz(-1.5547215) q[2];
sx q[2];
rz(-0.30259821) q[2];
rz(-0.042081632) q[3];
sx q[3];
rz(-0.29199568) q[3];
sx q[3];
rz(-0.80872768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06912032) q[0];
sx q[0];
rz(-2.0354249) q[0];
sx q[0];
rz(2.1492667) q[0];
rz(-2.8760257) q[1];
sx q[1];
rz(-2.3873603) q[1];
sx q[1];
rz(-3.1297562) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3385394) q[0];
sx q[0];
rz(-1.526471) q[0];
sx q[0];
rz(1.2903777) q[0];
rz(0.29143362) q[2];
sx q[2];
rz(-2.121392) q[2];
sx q[2];
rz(1.6432946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34764987) q[1];
sx q[1];
rz(-1.2041438) q[1];
sx q[1];
rz(-3.0701999) q[1];
rz(1.2740259) q[3];
sx q[3];
rz(-2.9057876) q[3];
sx q[3];
rz(0.91487003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.5806233) q[2];
sx q[2];
rz(-2.2443266) q[2];
sx q[2];
rz(2.238849) q[2];
rz(0.48815253) q[3];
sx q[3];
rz(-1.3267696) q[3];
sx q[3];
rz(-1.3443525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1667267) q[0];
sx q[0];
rz(-2.6868197) q[0];
sx q[0];
rz(1.7762666) q[0];
rz(2.4873554) q[1];
sx q[1];
rz(-0.74834329) q[1];
sx q[1];
rz(1.0985589) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64564686) q[0];
sx q[0];
rz(-0.80451999) q[0];
sx q[0];
rz(0.65973892) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6298348) q[2];
sx q[2];
rz(-1.8173886) q[2];
sx q[2];
rz(-2.5900019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2904086) q[1];
sx q[1];
rz(-2.1430204) q[1];
sx q[1];
rz(-1.6838724) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26596721) q[3];
sx q[3];
rz(-0.95925695) q[3];
sx q[3];
rz(-0.75436324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4174623) q[2];
sx q[2];
rz(-1.3516358) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2073479) q[0];
sx q[0];
rz(-1.8266015) q[0];
sx q[0];
rz(-2.9466002) q[0];
rz(-0.64942819) q[1];
sx q[1];
rz(-2.1959627) q[1];
sx q[1];
rz(1.6920998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3536384) q[0];
sx q[0];
rz(-2.1067841) q[0];
sx q[0];
rz(1.130442) q[0];
x q[1];
rz(-1.9282247) q[2];
sx q[2];
rz(-2.0783484) q[2];
sx q[2];
rz(-1.2032711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7342775) q[1];
sx q[1];
rz(-0.65562926) q[1];
sx q[1];
rz(-1.6796965) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4420801) q[3];
sx q[3];
rz(-1.8600924) q[3];
sx q[3];
rz(-0.6460748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0608757) q[2];
sx q[2];
rz(-2.1970811) q[2];
sx q[2];
rz(-0.057849217) q[2];
rz(0.079023376) q[3];
sx q[3];
rz(-1.9130324) q[3];
sx q[3];
rz(1.9144937) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3651315) q[0];
sx q[0];
rz(-1.4995898) q[0];
sx q[0];
rz(-0.32715964) q[0];
rz(-2.6311334) q[1];
sx q[1];
rz(-1.8233428) q[1];
sx q[1];
rz(-3.0697451) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851929) q[0];
sx q[0];
rz(-1.6677534) q[0];
sx q[0];
rz(0.18228874) q[0];
x q[1];
rz(-2.4113095) q[2];
sx q[2];
rz(-0.84368333) q[2];
sx q[2];
rz(-2.2471225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.209022) q[1];
sx q[1];
rz(-0.96197546) q[1];
sx q[1];
rz(0.8539416) q[1];
rz(-pi) q[2];
rz(1.9570146) q[3];
sx q[3];
rz(-0.90208331) q[3];
sx q[3];
rz(-2.6839601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0923126) q[2];
sx q[2];
rz(-1.4822373) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6698089) q[0];
sx q[0];
rz(-1.1265378) q[0];
sx q[0];
rz(-0.54927611) q[0];
rz(1.9688985) q[1];
sx q[1];
rz(-1.6212308) q[1];
sx q[1];
rz(1.8424013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0232721) q[0];
sx q[0];
rz(-2.9370312) q[0];
sx q[0];
rz(0.30318326) q[0];
x q[1];
rz(-1.9525398) q[2];
sx q[2];
rz(-0.92806584) q[2];
sx q[2];
rz(-2.1537154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0944984) q[1];
sx q[1];
rz(-0.21378042) q[1];
sx q[1];
rz(0.70014145) q[1];
rz(-0.48051832) q[3];
sx q[3];
rz(-2.5439265) q[3];
sx q[3];
rz(-1.0880119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4486763) q[2];
sx q[2];
rz(-2.0557949) q[2];
sx q[2];
rz(-0.14796251) q[2];
rz(2.1678917) q[3];
sx q[3];
rz(-2.2805043) q[3];
sx q[3];
rz(-0.93419832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14514087) q[0];
sx q[0];
rz(-1.3057764) q[0];
sx q[0];
rz(1.8291352) q[0];
rz(1.3120069) q[1];
sx q[1];
rz(-0.94999718) q[1];
sx q[1];
rz(-2.2012262) q[1];
rz(-1.1205705) q[2];
sx q[2];
rz(-2.8489589) q[2];
sx q[2];
rz(2.6363303) q[2];
rz(1.3668904) q[3];
sx q[3];
rz(-0.6253607) q[3];
sx q[3];
rz(-1.6224773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
