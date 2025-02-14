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
rz(1.4580392) q[0];
sx q[0];
rz(-0.91705051) q[0];
sx q[0];
rz(0.24669692) q[0];
rz(-2.2973581) q[1];
sx q[1];
rz(-1.8480453) q[1];
sx q[1];
rz(-1.3026613) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5555252) q[0];
sx q[0];
rz(-1.925615) q[0];
sx q[0];
rz(0.3094425) q[0];
rz(-2.1470492) q[2];
sx q[2];
rz(-0.79205293) q[2];
sx q[2];
rz(-2.920237) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9907551) q[1];
sx q[1];
rz(-2.1937943) q[1];
sx q[1];
rz(1.4332818) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0329923) q[3];
sx q[3];
rz(-1.7157156) q[3];
sx q[3];
rz(2.8371881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.777433) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(-1.1945266) q[2];
rz(1.2398237) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(-1.4136081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164923) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(-0.66705739) q[0];
rz(-0.27101135) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(1.0557231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8264814) q[0];
sx q[0];
rz(-1.8391063) q[0];
sx q[0];
rz(0.92197355) q[0];
rz(-3.1201601) q[2];
sx q[2];
rz(-1.9551245) q[2];
sx q[2];
rz(0.0057980428) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.70433863) q[1];
sx q[1];
rz(-1.9074215) q[1];
sx q[1];
rz(3.047154) q[1];
rz(-pi) q[2];
rz(-0.74774489) q[3];
sx q[3];
rz(-0.90255957) q[3];
sx q[3];
rz(-1.8577607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6229652) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(0.61678994) q[2];
rz(0.24584298) q[3];
sx q[3];
rz(-1.5522233) q[3];
sx q[3];
rz(-1.8776548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.371405) q[0];
sx q[0];
rz(-0.73391947) q[0];
sx q[0];
rz(-0.1804633) q[0];
rz(2.6626185) q[1];
sx q[1];
rz(-2.6909747) q[1];
sx q[1];
rz(-1.825038) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0947803) q[0];
sx q[0];
rz(-0.98623305) q[0];
sx q[0];
rz(0.78740904) q[0];
x q[1];
rz(-2.0094078) q[2];
sx q[2];
rz(-0.45044611) q[2];
sx q[2];
rz(0.23368719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7261947) q[1];
sx q[1];
rz(-1.2737927) q[1];
sx q[1];
rz(2.4261977) q[1];
rz(-pi) q[2];
rz(2.4495921) q[3];
sx q[3];
rz(-2.3131136) q[3];
sx q[3];
rz(0.52317515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58671826) q[2];
sx q[2];
rz(-2.6719533) q[2];
sx q[2];
rz(-0.47373104) q[2];
rz(-2.6321865) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(1.1562851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8944775) q[0];
sx q[0];
rz(-3.101427) q[0];
sx q[0];
rz(-1.0152869) q[0];
rz(-2.9292551) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(2.7755348) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057945874) q[0];
sx q[0];
rz(-2.2444631) q[0];
sx q[0];
rz(0.64274733) q[0];
rz(-pi) q[1];
rz(-2.9061142) q[2];
sx q[2];
rz(-1.7120541) q[2];
sx q[2];
rz(1.3154674) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8180994) q[1];
sx q[1];
rz(-0.93902367) q[1];
sx q[1];
rz(-1.9629823) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96694209) q[3];
sx q[3];
rz(-1.8181385) q[3];
sx q[3];
rz(0.027146904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66146835) q[2];
sx q[2];
rz(-1.3084359) q[2];
sx q[2];
rz(-0.69501957) q[2];
rz(-0.85092893) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(-1.883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6224391) q[0];
sx q[0];
rz(-0.16534403) q[0];
sx q[0];
rz(-2.8261321) q[0];
rz(2.8568965) q[1];
sx q[1];
rz(-1.618914) q[1];
sx q[1];
rz(-0.42627898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5770059) q[0];
sx q[0];
rz(-0.88250151) q[0];
sx q[0];
rz(0.37969638) q[0];
rz(-pi) q[1];
rz(-2.2441909) q[2];
sx q[2];
rz(-0.36172141) q[2];
sx q[2];
rz(0.42332403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46968383) q[1];
sx q[1];
rz(-1.4844746) q[1];
sx q[1];
rz(1.1344019) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0067389) q[3];
sx q[3];
rz(-2.5814179) q[3];
sx q[3];
rz(-0.87040802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2686501) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(3.1128913) q[2];
rz(-0.21102333) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(-2.4215462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215304) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(2.9826214) q[0];
rz(0.65525118) q[1];
sx q[1];
rz(-2.7924004) q[1];
sx q[1];
rz(-1.6045301) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6626016) q[0];
sx q[0];
rz(-2.3259386) q[0];
sx q[0];
rz(-2.2032477) q[0];
rz(-pi) q[1];
rz(-0.75266204) q[2];
sx q[2];
rz(-0.48303451) q[2];
sx q[2];
rz(2.8217884) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3517396) q[1];
sx q[1];
rz(-2.5765214) q[1];
sx q[1];
rz(-1.5365824) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37015583) q[3];
sx q[3];
rz(-1.2803382) q[3];
sx q[3];
rz(-0.53122093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0752461) q[2];
sx q[2];
rz(-1.9074351) q[2];
sx q[2];
rz(0.1417024) q[2];
rz(3.1250478) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(2.580548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4905106) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(0.91019994) q[0];
rz(0.74288145) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(1.9076617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54523477) q[0];
sx q[0];
rz(-1.4995793) q[0];
sx q[0];
rz(-0.13990732) q[0];
x q[1];
rz(2.874159) q[2];
sx q[2];
rz(-1.1743011) q[2];
sx q[2];
rz(-0.948179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9561718) q[1];
sx q[1];
rz(-2.2347576) q[1];
sx q[1];
rz(-1.0770304) q[1];
x q[2];
rz(2.5508826) q[3];
sx q[3];
rz(-1.9001257) q[3];
sx q[3];
rz(-1.2423837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10741216) q[2];
sx q[2];
rz(-1.567652) q[2];
sx q[2];
rz(0.91599715) q[2];
rz(0.25137526) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(0.34534064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4995572) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(3.1009951) q[0];
rz(-1.1740855) q[1];
sx q[1];
rz(-0.65316713) q[1];
sx q[1];
rz(-2.0814799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030329372) q[0];
sx q[0];
rz(-1.666774) q[0];
sx q[0];
rz(2.7594271) q[0];
x q[1];
rz(1.4547075) q[2];
sx q[2];
rz(-1.5512067) q[2];
sx q[2];
rz(1.3649114) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6190784) q[1];
sx q[1];
rz(-2.1685617) q[1];
sx q[1];
rz(1.9512216) q[1];
x q[2];
rz(-2.4828217) q[3];
sx q[3];
rz(-2.3537354) q[3];
sx q[3];
rz(3.0751021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4982831) q[2];
sx q[2];
rz(-2.0556367) q[2];
sx q[2];
rz(2.3504284) q[2];
rz(1.5353954) q[3];
sx q[3];
rz(-0.94208661) q[3];
sx q[3];
rz(0.78996381) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20752792) q[0];
sx q[0];
rz(-2.9495033) q[0];
sx q[0];
rz(2.9839363) q[0];
rz(-0.12570307) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(-0.30473614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5253741) q[0];
sx q[0];
rz(-2.8044639) q[0];
sx q[0];
rz(1.3167156) q[0];
x q[1];
rz(-1.1515806) q[2];
sx q[2];
rz(-2.1682924) q[2];
sx q[2];
rz(-2.1066163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8915776) q[1];
sx q[1];
rz(-3.0072837) q[1];
sx q[1];
rz(3.1242643) q[1];
rz(2.7823506) q[3];
sx q[3];
rz(-0.83165681) q[3];
sx q[3];
rz(-3.0225282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8704845) q[2];
sx q[2];
rz(-1.5616337) q[2];
sx q[2];
rz(-1.6893207) q[2];
rz(-2.3952386) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-0.11317429) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8993503) q[0];
sx q[0];
rz(-0.32651383) q[0];
sx q[0];
rz(0.53264701) q[0];
rz(0.38231725) q[1];
sx q[1];
rz(-2.2492354) q[1];
sx q[1];
rz(-0.16758448) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78412752) q[0];
sx q[0];
rz(-1.3793769) q[0];
sx q[0];
rz(-3.0166059) q[0];
rz(-pi) q[1];
rz(-0.26248502) q[2];
sx q[2];
rz(-1.0392611) q[2];
sx q[2];
rz(2.558311) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6217612) q[1];
sx q[1];
rz(-1.863136) q[1];
sx q[1];
rz(1.2755434) q[1];
rz(-pi) q[2];
rz(-0.83491171) q[3];
sx q[3];
rz(-1.8855699) q[3];
sx q[3];
rz(1.9251458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8297537) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(-0.46693841) q[2];
rz(-1.1030819) q[3];
sx q[3];
rz(-1.1934049) q[3];
sx q[3];
rz(-1.232049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5760096) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(-0.78454984) q[1];
sx q[1];
rz(-0.34930925) q[1];
sx q[1];
rz(0.92225155) q[1];
rz(2.2837737) q[2];
sx q[2];
rz(-1.0150649) q[2];
sx q[2];
rz(-0.82761717) q[2];
rz(1.2123232) q[3];
sx q[3];
rz(-1.3176821) q[3];
sx q[3];
rz(-1.1715922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
