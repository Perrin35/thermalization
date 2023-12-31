OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(2.1587125) q[0];
sx q[0];
rz(8.4150253) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(1.8523822) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56611094) q[0];
sx q[0];
rz(-2.5876343) q[0];
sx q[0];
rz(-2.1411116) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5640366) q[2];
sx q[2];
rz(-1.5951556) q[2];
sx q[2];
rz(-1.3855795) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7154555) q[1];
sx q[1];
rz(-1.8278367) q[1];
sx q[1];
rz(-2.9261158) q[1];
x q[2];
rz(-0.07006499) q[3];
sx q[3];
rz(-1.3032039) q[3];
sx q[3];
rz(-2.23578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5499128) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(-3.0060449) q[2];
rz(0.1401976) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(-0.26023284) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(1.82812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03598729) q[0];
sx q[0];
rz(-0.30788883) q[0];
sx q[0];
rz(2.3769828) q[0];
x q[1];
rz(2.1140695) q[2];
sx q[2];
rz(-2.5894905) q[2];
sx q[2];
rz(-3.0039624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62994176) q[1];
sx q[1];
rz(-0.18254666) q[1];
sx q[1];
rz(2.6444142) q[1];
rz(-pi) q[2];
rz(0.95150681) q[3];
sx q[3];
rz(-3.0118239) q[3];
sx q[3];
rz(1.4459923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69497067) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(2.4296956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057782877) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-0.69586786) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(2.9755039) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9244708) q[0];
sx q[0];
rz(-2.8170601) q[0];
sx q[0];
rz(1.0990418) q[0];
rz(2.5327352) q[2];
sx q[2];
rz(-1.1470084) q[2];
sx q[2];
rz(-0.54178967) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3661763) q[1];
sx q[1];
rz(-2.011538) q[1];
sx q[1];
rz(-2.612545) q[1];
rz(-pi) q[2];
rz(1.2766248) q[3];
sx q[3];
rz(-2.4066381) q[3];
sx q[3];
rz(0.24925772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(0.44149533) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7168717) q[0];
sx q[0];
rz(-1.4972591) q[0];
sx q[0];
rz(1.5872692) q[0];
x q[1];
rz(1.5261126) q[2];
sx q[2];
rz(-1.2031021) q[2];
sx q[2];
rz(1.1553264) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2312154) q[1];
sx q[1];
rz(-2.7312355) q[1];
sx q[1];
rz(-1.7675722) q[1];
rz(-pi) q[2];
rz(2.1500548) q[3];
sx q[3];
rz(-2.8354128) q[3];
sx q[3];
rz(-1.7508208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(1.305497) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(1.501804) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.775979) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-0.10770527) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90785039) q[0];
sx q[0];
rz(-1.478501) q[0];
sx q[0];
rz(1.924563) q[0];
rz(-pi) q[1];
rz(-0.55107112) q[2];
sx q[2];
rz(-2.1398009) q[2];
sx q[2];
rz(2.0713194) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2376033) q[1];
sx q[1];
rz(-0.9346107) q[1];
sx q[1];
rz(-2.643305) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7961736) q[3];
sx q[3];
rz(-2.1072227) q[3];
sx q[3];
rz(0.1766583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0832131) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(1.1389114) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(-1.0169792) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6769584) q[0];
sx q[0];
rz(-0.58566739) q[0];
sx q[0];
rz(-2.2023489) q[0];
rz(0.20118841) q[2];
sx q[2];
rz(-1.5474461) q[2];
sx q[2];
rz(0.77229283) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.444126) q[1];
sx q[1];
rz(-1.0870458) q[1];
sx q[1];
rz(-2.6078348) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0820886) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(0.90261501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(0.44657782) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25062659) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.9218146) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0735002) q[0];
sx q[0];
rz(-1.6627899) q[0];
sx q[0];
rz(-1.5324701) q[0];
rz(-pi) q[1];
rz(2.2336823) q[2];
sx q[2];
rz(-2.6320576) q[2];
sx q[2];
rz(0.31928911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4525758) q[1];
sx q[1];
rz(-2.2170482) q[1];
sx q[1];
rz(-2.8313682) q[1];
x q[2];
rz(1.9795498) q[3];
sx q[3];
rz(-1.4882653) q[3];
sx q[3];
rz(-2.2857034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(-3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(-2.3855456) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379631) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(-1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(-2.7323639) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17657875) q[0];
sx q[0];
rz(-1.8518378) q[0];
sx q[0];
rz(0.9126419) q[0];
rz(-pi) q[1];
rz(-2.6312469) q[2];
sx q[2];
rz(-0.58594698) q[2];
sx q[2];
rz(2.3454391) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.65539) q[1];
sx q[1];
rz(-1.4185925) q[1];
sx q[1];
rz(-0.57032077) q[1];
rz(2.6243022) q[3];
sx q[3];
rz(-1.358277) q[3];
sx q[3];
rz(-1.1503435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(0.64152843) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(2.0857281) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(-1.508629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2471837) q[0];
sx q[0];
rz(-1.7315454) q[0];
sx q[0];
rz(-2.2649691) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0502003) q[2];
sx q[2];
rz(-0.95801991) q[2];
sx q[2];
rz(-2.6932655) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0548045) q[1];
sx q[1];
rz(-2.0874223) q[1];
sx q[1];
rz(-0.36887849) q[1];
x q[2];
rz(0.86767254) q[3];
sx q[3];
rz(-1.5276067) q[3];
sx q[3];
rz(-1.6004576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325539) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(-2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(-1.7601097) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59795838) q[0];
sx q[0];
rz(-1.6126452) q[0];
sx q[0];
rz(-2.6393487) q[0];
x q[1];
rz(-2.9276766) q[2];
sx q[2];
rz(-2.6579318) q[2];
sx q[2];
rz(0.96825251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46422058) q[1];
sx q[1];
rz(-0.83253011) q[1];
sx q[1];
rz(-0.99454753) q[1];
rz(-pi) q[2];
rz(2.8548106) q[3];
sx q[3];
rz(-2.2638595) q[3];
sx q[3];
rz(1.9758488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.4505724) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.272841) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(-1.8252774) q[2];
sx q[2];
rz(-1.5263867) q[2];
sx q[2];
rz(2.9771794) q[2];
rz(0.5100526) q[3];
sx q[3];
rz(-1.4997713) q[3];
sx q[3];
rz(-0.98239519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
