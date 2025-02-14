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
rz(3.0348294) q[0];
sx q[0];
rz(-1.3314629) q[0];
sx q[0];
rz(-1.7911628) q[0];
rz(-2.8332233) q[1];
sx q[1];
rz(-2.8928533) q[1];
sx q[1];
rz(-2.1623936) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0672507) q[0];
sx q[0];
rz(-2.3789211) q[0];
sx q[0];
rz(-3.0347155) q[0];
rz(-pi) q[1];
rz(-1.1264231) q[2];
sx q[2];
rz(-2.4326486) q[2];
sx q[2];
rz(-2.5042748) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1147656) q[1];
sx q[1];
rz(-1.9239167) q[1];
sx q[1];
rz(1.4985282) q[1];
x q[2];
rz(0.42433831) q[3];
sx q[3];
rz(-2.7380214) q[3];
sx q[3];
rz(-0.85116019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.2110445) q[2];
sx q[2];
rz(-0.63300526) q[2];
rz(2.4999319) q[3];
sx q[3];
rz(-2.6810665) q[3];
sx q[3];
rz(0.7707001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68384701) q[0];
sx q[0];
rz(-0.96413833) q[0];
sx q[0];
rz(-2.316851) q[0];
rz(1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(2.8214084) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31486785) q[0];
sx q[0];
rz(-1.1629675) q[0];
sx q[0];
rz(2.639222) q[0];
rz(-1.0616706) q[2];
sx q[2];
rz(-1.7560282) q[2];
sx q[2];
rz(2.9502466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0517438) q[1];
sx q[1];
rz(-1.4898572) q[1];
sx q[1];
rz(-1.5481987) q[1];
x q[2];
rz(-1.792644) q[3];
sx q[3];
rz(-1.5127521) q[3];
sx q[3];
rz(-2.0741418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0188633) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(1.2790206) q[2];
rz(3.0458798) q[3];
sx q[3];
rz(-1.0749823) q[3];
sx q[3];
rz(1.3191282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9886446) q[0];
sx q[0];
rz(-2.4929292) q[0];
sx q[0];
rz(-1.0308107) q[0];
rz(-0.55054322) q[1];
sx q[1];
rz(-2.2830453) q[1];
sx q[1];
rz(1.2446838) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1079946) q[0];
sx q[0];
rz(-1.4116086) q[0];
sx q[0];
rz(1.6248996) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95374505) q[2];
sx q[2];
rz(-1.707952) q[2];
sx q[2];
rz(2.0037665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.479949) q[1];
sx q[1];
rz(-2.2315761) q[1];
sx q[1];
rz(-0.26008545) q[1];
rz(-pi) q[2];
rz(-2.3233285) q[3];
sx q[3];
rz(-1.8570559) q[3];
sx q[3];
rz(0.32901007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9024258) q[2];
sx q[2];
rz(-0.5639762) q[2];
sx q[2];
rz(-1.172056) q[2];
rz(0.82550448) q[3];
sx q[3];
rz(-1.2647311) q[3];
sx q[3];
rz(0.66506213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.091752) q[0];
sx q[0];
rz(-1.4565775) q[0];
sx q[0];
rz(-2.7276584) q[0];
rz(2.6992758) q[1];
sx q[1];
rz(-1.0041753) q[1];
sx q[1];
rz(0.54534674) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2065679) q[0];
sx q[0];
rz(-0.27982084) q[0];
sx q[0];
rz(2.9601475) q[0];
x q[1];
rz(-2.3661883) q[2];
sx q[2];
rz(-1.0808965) q[2];
sx q[2];
rz(2.6012914) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87121039) q[1];
sx q[1];
rz(-2.5426513) q[1];
sx q[1];
rz(-1.0652871) q[1];
rz(-pi) q[2];
rz(2.993587) q[3];
sx q[3];
rz(-2.0750351) q[3];
sx q[3];
rz(-1.3324236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.44117323) q[2];
sx q[2];
rz(-1.2607231) q[2];
sx q[2];
rz(-2.9803661) q[2];
rz(1.5876596) q[3];
sx q[3];
rz(-2.5416538) q[3];
sx q[3];
rz(0.59066311) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0597543) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(-2.4181714) q[0];
rz(1.8266034) q[1];
sx q[1];
rz(-1.864121) q[1];
sx q[1];
rz(-2.1037197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4525131) q[0];
sx q[0];
rz(-2.291579) q[0];
sx q[0];
rz(1.4656653) q[0];
rz(-pi) q[1];
rz(2.5332301) q[2];
sx q[2];
rz(-2.3181097) q[2];
sx q[2];
rz(2.0840816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.65601634) q[1];
sx q[1];
rz(-0.52158251) q[1];
sx q[1];
rz(0.55934577) q[1];
rz(-pi) q[2];
rz(2.5955022) q[3];
sx q[3];
rz(-1.6973064) q[3];
sx q[3];
rz(-1.5777335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7448685) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(-0.57207668) q[2];
rz(0.62503302) q[3];
sx q[3];
rz(-2.1376938) q[3];
sx q[3];
rz(1.9151789) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52142757) q[0];
sx q[0];
rz(-0.033997424) q[0];
sx q[0];
rz(0.52325621) q[0];
rz(-1.322586) q[1];
sx q[1];
rz(-1.6736284) q[1];
sx q[1];
rz(2.5214213) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5145743) q[0];
sx q[0];
rz(-0.73643273) q[0];
sx q[0];
rz(-2.4798142) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6441395) q[2];
sx q[2];
rz(-1.1032602) q[2];
sx q[2];
rz(-0.95638025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7249696) q[1];
sx q[1];
rz(-2.1955804) q[1];
sx q[1];
rz(-0.10819541) q[1];
x q[2];
rz(-1.3955388) q[3];
sx q[3];
rz(-1.3564928) q[3];
sx q[3];
rz(-2.8068723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4124477) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(1.283851) q[2];
rz(2.1733952) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(-1.5549972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9667483) q[0];
sx q[0];
rz(-1.1356249) q[0];
sx q[0];
rz(-1.9298166) q[0];
rz(2.1719596) q[1];
sx q[1];
rz(-1.8200834) q[1];
sx q[1];
rz(1.8671573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2454555) q[0];
sx q[0];
rz(-1.7217759) q[0];
sx q[0];
rz(0.87489382) q[0];
rz(-pi) q[1];
rz(0.095000141) q[2];
sx q[2];
rz(-2.3534905) q[2];
sx q[2];
rz(-2.1301477) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.798118) q[1];
sx q[1];
rz(-0.96303446) q[1];
sx q[1];
rz(-0.37109427) q[1];
rz(-pi) q[2];
rz(1.2479374) q[3];
sx q[3];
rz(-2.0104694) q[3];
sx q[3];
rz(0.46931258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6756639) q[2];
sx q[2];
rz(-1.6589386) q[2];
sx q[2];
rz(-2.8181804) q[2];
rz(-0.0023500738) q[3];
sx q[3];
rz(-1.4582783) q[3];
sx q[3];
rz(2.5200444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49331409) q[0];
sx q[0];
rz(-1.4819773) q[0];
sx q[0];
rz(1.3508654) q[0];
rz(-1.3510652) q[1];
sx q[1];
rz(-1.8565145) q[1];
sx q[1];
rz(-1.2877119) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5918811) q[0];
sx q[0];
rz(-1.2882107) q[0];
sx q[0];
rz(-3.0192399) q[0];
rz(-pi) q[1];
rz(0.28952629) q[2];
sx q[2];
rz(-2.2643201) q[2];
sx q[2];
rz(-2.5159581) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2877174) q[1];
sx q[1];
rz(-0.36810222) q[1];
sx q[1];
rz(-0.48990695) q[1];
rz(-pi) q[2];
rz(-1.5519556) q[3];
sx q[3];
rz(-1.9196438) q[3];
sx q[3];
rz(-0.89295372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52435851) q[2];
sx q[2];
rz(-1.1523767) q[2];
sx q[2];
rz(2.1481245) q[2];
rz(-1.0348381) q[3];
sx q[3];
rz(-0.60641685) q[3];
sx q[3];
rz(2.9724227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33191037) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(1.6140953) q[0];
rz(2.2612259) q[1];
sx q[1];
rz(-0.4207193) q[1];
sx q[1];
rz(0.31164935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0535705) q[0];
sx q[0];
rz(-0.35614518) q[0];
sx q[0];
rz(0.97266622) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8021605) q[2];
sx q[2];
rz(-1.5867434) q[2];
sx q[2];
rz(1.1167662) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91395177) q[1];
sx q[1];
rz(-0.80004179) q[1];
sx q[1];
rz(2.6197919) q[1];
rz(2.8086189) q[3];
sx q[3];
rz(-0.94374386) q[3];
sx q[3];
rz(2.0505333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6722022) q[2];
sx q[2];
rz(-2.2180874) q[2];
sx q[2];
rz(-1.3908609) q[2];
rz(-1.5270799) q[3];
sx q[3];
rz(-2.0938087) q[3];
sx q[3];
rz(1.0745777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5020318) q[0];
sx q[0];
rz(-1.1447516) q[0];
sx q[0];
rz(2.5290053) q[0];
rz(-2.1514429) q[1];
sx q[1];
rz(-1.1724816) q[1];
sx q[1];
rz(1.6506857) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6698786) q[0];
sx q[0];
rz(-1.9900787) q[0];
sx q[0];
rz(-0.74434481) q[0];
rz(-1.2464929) q[2];
sx q[2];
rz(-0.68585472) q[2];
sx q[2];
rz(1.7795479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1742434) q[1];
sx q[1];
rz(-0.16160204) q[1];
sx q[1];
rz(1.5320278) q[1];
x q[2];
rz(2.4435639) q[3];
sx q[3];
rz(-1.7637611) q[3];
sx q[3];
rz(-0.86856996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7154197) q[2];
sx q[2];
rz(-0.82948589) q[2];
sx q[2];
rz(-0.034991525) q[2];
rz(-1.70111) q[3];
sx q[3];
rz(-0.8927497) q[3];
sx q[3];
rz(0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424425) q[0];
sx q[0];
rz(-2.2461666) q[0];
sx q[0];
rz(2.9641892) q[0];
rz(1.9433446) q[1];
sx q[1];
rz(-1.6234963) q[1];
sx q[1];
rz(-2.1877098) q[1];
rz(-2.162355) q[2];
sx q[2];
rz(-2.1605347) q[2];
sx q[2];
rz(0.84183358) q[2];
rz(2.421834) q[3];
sx q[3];
rz(-1.2133141) q[3];
sx q[3];
rz(2.6837466) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
