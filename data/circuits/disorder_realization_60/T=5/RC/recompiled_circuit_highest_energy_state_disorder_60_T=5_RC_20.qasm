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
rz(-1.6835535) q[0];
sx q[0];
rz(-2.2245421) q[0];
sx q[0];
rz(2.8948957) q[0];
rz(-2.2973581) q[1];
sx q[1];
rz(-1.8480453) q[1];
sx q[1];
rz(1.8389314) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8116353) q[0];
sx q[0];
rz(-0.46644743) q[0];
sx q[0];
rz(-2.2587725) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1470492) q[2];
sx q[2];
rz(-0.79205293) q[2];
sx q[2];
rz(-2.920237) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9907551) q[1];
sx q[1];
rz(-2.1937943) q[1];
sx q[1];
rz(1.7083108) q[1];
rz(-pi) q[2];
rz(-1.4250303) q[3];
sx q[3];
rz(-1.4633388) q[3];
sx q[3];
rz(-1.859457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.777433) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(1.1945266) q[2];
rz(-1.901769) q[3];
sx q[3];
rz(-1.4907962) q[3];
sx q[3];
rz(1.4136081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0251004) q[0];
sx q[0];
rz(-1.1273552) q[0];
sx q[0];
rz(0.66705739) q[0];
rz(2.8705813) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(-2.0858696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0572898) q[0];
sx q[0];
rz(-0.94883666) q[0];
sx q[0];
rz(0.33226407) q[0];
x q[1];
rz(1.6237359) q[2];
sx q[2];
rz(-2.756697) q[2];
sx q[2];
rz(-0.062907779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.437254) q[1];
sx q[1];
rz(-1.9074215) q[1];
sx q[1];
rz(3.047154) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74774489) q[3];
sx q[3];
rz(-0.90255957) q[3];
sx q[3];
rz(-1.2838319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6229652) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(-0.61678994) q[2];
rz(-0.24584298) q[3];
sx q[3];
rz(-1.5522233) q[3];
sx q[3];
rz(1.8776548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7701876) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(2.9611294) q[0];
rz(0.47897419) q[1];
sx q[1];
rz(-0.45061794) q[1];
sx q[1];
rz(-1.825038) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265357) q[0];
sx q[0];
rz(-2.2000072) q[0];
sx q[0];
rz(2.3903484) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9834824) q[2];
sx q[2];
rz(-1.3848334) q[2];
sx q[2];
rz(0.93753147) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2352202) q[1];
sx q[1];
rz(-0.89284578) q[1];
sx q[1];
rz(1.1855907) q[1];
rz(-pi) q[2];
rz(0.69200055) q[3];
sx q[3];
rz(-0.82847906) q[3];
sx q[3];
rz(0.52317515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.58671826) q[2];
sx q[2];
rz(-0.46963936) q[2];
sx q[2];
rz(-0.47373104) q[2];
rz(0.50940618) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(-1.9853076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944775) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(-1.0152869) q[0];
rz(-2.9292551) q[1];
sx q[1];
rz(-1.5539853) q[1];
sx q[1];
rz(0.36605787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81824077) q[0];
sx q[0];
rz(-0.89491544) q[0];
sx q[0];
rz(0.9267207) q[0];
x q[1];
rz(2.9061142) q[2];
sx q[2];
rz(-1.4295385) q[2];
sx q[2];
rz(-1.8261253) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4868813) q[1];
sx q[1];
rz(-1.2572479) q[1];
sx q[1];
rz(-2.4717719) q[1];
x q[2];
rz(0.96694209) q[3];
sx q[3];
rz(-1.8181385) q[3];
sx q[3];
rz(0.027146904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66146835) q[2];
sx q[2];
rz(-1.3084359) q[2];
sx q[2];
rz(-0.69501957) q[2];
rz(2.2906637) q[3];
sx q[3];
rz(-1.3635819) q[3];
sx q[3];
rz(-1.2581576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6224391) q[0];
sx q[0];
rz(-0.16534403) q[0];
sx q[0];
rz(0.31546053) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(-2.7153137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56458679) q[0];
sx q[0];
rz(-2.2590911) q[0];
sx q[0];
rz(-2.7618963) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23172863) q[2];
sx q[2];
rz(-1.2905057) q[2];
sx q[2];
rz(-2.0120399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0002928) q[1];
sx q[1];
rz(-2.0054549) q[1];
sx q[1];
rz(3.0463957) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0067389) q[3];
sx q[3];
rz(-0.56017471) q[3];
sx q[3];
rz(0.87040802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8729426) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(0.028701393) q[2];
rz(-2.9305693) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(-0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215304) q[0];
sx q[0];
rz(-1.2380607) q[0];
sx q[0];
rz(-2.9826214) q[0];
rz(0.65525118) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(1.6045301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6626016) q[0];
sx q[0];
rz(-2.3259386) q[0];
sx q[0];
rz(0.93834491) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9150429) q[2];
sx q[2];
rz(-1.2249399) q[2];
sx q[2];
rz(2.648022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83035806) q[1];
sx q[1];
rz(-1.0060961) q[1];
sx q[1];
rz(3.1199074) q[1];
rz(-pi) q[2];
rz(-0.69055478) q[3];
sx q[3];
rz(-0.46633807) q[3];
sx q[3];
rz(0.40377221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0752461) q[2];
sx q[2];
rz(-1.9074351) q[2];
sx q[2];
rz(0.1417024) q[2];
rz(-0.016544841) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(2.580548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510821) q[0];
sx q[0];
rz(-0.53618479) q[0];
sx q[0];
rz(0.91019994) q[0];
rz(2.3987112) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(-1.9076617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5963579) q[0];
sx q[0];
rz(-1.6420134) q[0];
sx q[0];
rz(-0.13990732) q[0];
rz(-pi) q[1];
rz(-1.0077613) q[2];
sx q[2];
rz(-0.47426155) q[2];
sx q[2];
rz(-0.33111085) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.90335315) q[1];
sx q[1];
rz(-0.80447865) q[1];
sx q[1];
rz(0.54460214) q[1];
x q[2];
rz(-0.55039946) q[3];
sx q[3];
rz(-0.66662753) q[3];
sx q[3];
rz(2.3638099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0341805) q[2];
sx q[2];
rz(-1.567652) q[2];
sx q[2];
rz(2.2255955) q[2];
rz(2.8902174) q[3];
sx q[3];
rz(-2.1706457) q[3];
sx q[3];
rz(-2.796252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6420355) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(-0.04059759) q[0];
rz(-1.1740855) q[1];
sx q[1];
rz(-0.65316713) q[1];
sx q[1];
rz(1.0601128) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3670676) q[0];
sx q[0];
rz(-2.7481348) q[0];
sx q[0];
rz(2.8889546) q[0];
x q[1];
rz(0.019722299) q[2];
sx q[2];
rz(-1.6868627) q[2];
sx q[2];
rz(-2.937992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0012944) q[1];
sx q[1];
rz(-2.4457275) q[1];
sx q[1];
rz(-2.6423179) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4701263) q[3];
sx q[3];
rz(-2.0196303) q[3];
sx q[3];
rz(-2.0041114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4982831) q[2];
sx q[2];
rz(-1.085956) q[2];
sx q[2];
rz(0.79116428) q[2];
rz(1.6061973) q[3];
sx q[3];
rz(-2.199506) q[3];
sx q[3];
rz(-2.3516288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20752792) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(2.9839363) q[0];
rz(3.0158896) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(-0.30473614) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9466772) q[0];
sx q[0];
rz(-1.4875571) q[0];
sx q[0];
rz(-1.2437245) q[0];
rz(-pi) q[1];
rz(-2.6025099) q[2];
sx q[2];
rz(-2.4266908) q[2];
sx q[2];
rz(1.4366921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23252919) q[1];
sx q[1];
rz(-1.705085) q[1];
sx q[1];
rz(1.5731377) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35924201) q[3];
sx q[3];
rz(-2.3099358) q[3];
sx q[3];
rz(0.11906448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8704845) q[2];
sx q[2];
rz(-1.5799589) q[2];
sx q[2];
rz(1.6893207) q[2];
rz(-2.3952386) q[3];
sx q[3];
rz(-1.690026) q[3];
sx q[3];
rz(-3.0284184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24224237) q[0];
sx q[0];
rz(-2.8150788) q[0];
sx q[0];
rz(-0.53264701) q[0];
rz(0.38231725) q[1];
sx q[1];
rz(-0.89235726) q[1];
sx q[1];
rz(0.16758448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9411128) q[0];
sx q[0];
rz(-0.22819209) q[0];
sx q[0];
rz(0.99916332) q[0];
rz(-pi) q[1];
rz(-1.0239059) q[2];
sx q[2];
rz(-1.3452072) q[2];
sx q[2];
rz(1.1228665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1780562) q[1];
sx q[1];
rz(-1.85317) q[1];
sx q[1];
rz(2.8368188) q[1];
rz(-pi) q[2];
rz(2.3066809) q[3];
sx q[3];
rz(-1.2560227) q[3];
sx q[3];
rz(1.2164468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.31183895) q[2];
sx q[2];
rz(-1.6646174) q[2];
sx q[2];
rz(2.6746542) q[2];
rz(1.1030819) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(1.9095437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5655831) q[0];
sx q[0];
rz(-0.82721114) q[0];
sx q[0];
rz(0.4785434) q[0];
rz(2.3570428) q[1];
sx q[1];
rz(-0.34930925) q[1];
sx q[1];
rz(0.92225155) q[1];
rz(2.4541773) q[2];
sx q[2];
rz(-0.981642) q[2];
sx q[2];
rz(-2.8264075) q[2];
rz(-1.2123232) q[3];
sx q[3];
rz(-1.8239106) q[3];
sx q[3];
rz(1.9700005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
