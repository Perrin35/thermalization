OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4607853) q[0];
sx q[0];
rz(-2.1587125) q[0];
sx q[0];
rz(2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(1.8523822) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56611094) q[0];
sx q[0];
rz(-0.55395836) q[0];
sx q[0];
rz(-1.0004811) q[0];
rz(-2.870954) q[2];
sx q[2];
rz(-0.025279609) q[2];
sx q[2];
rz(1.6563005) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28469052) q[1];
sx q[1];
rz(-2.8077217) q[1];
sx q[1];
rz(-2.2536709) q[1];
x q[2];
rz(1.8390158) q[3];
sx q[3];
rz(-1.5032288) q[3];
sx q[3];
rz(-0.64642954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(0.1401976) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(-2.360789) q[0];
rz(-0.26023284) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.82812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03598729) q[0];
sx q[0];
rz(-2.8337038) q[0];
sx q[0];
rz(2.3769828) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30828373) q[2];
sx q[2];
rz(-1.1051902) q[2];
sx q[2];
rz(-0.75454933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7105512) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(2.980742) q[1];
rz(-pi) q[2];
rz(1.6766657) q[3];
sx q[3];
rz(-1.6459811) q[3];
sx q[3];
rz(0.74010805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(1.4387087) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(0.71189705) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057782877) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(0.16608873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9244708) q[0];
sx q[0];
rz(-2.8170601) q[0];
sx q[0];
rz(-2.0425509) q[0];
rz(-pi) q[1];
rz(-0.66785779) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(-2.6454676) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6925466) q[1];
sx q[1];
rz(-2.0448301) q[1];
sx q[1];
rz(-1.0707335) q[1];
rz(-pi) q[2];
rz(0.8576287) q[3];
sx q[3];
rz(-1.3751251) q[3];
sx q[3];
rz(1.5426202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(1.123547) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0825901) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(0.77484432) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42472096) q[0];
sx q[0];
rz(-1.4972591) q[0];
sx q[0];
rz(1.5872692) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7735633) q[2];
sx q[2];
rz(-1.5291011) q[2];
sx q[2];
rz(-0.3993984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0171417) q[1];
sx q[1];
rz(-1.1688197) q[1];
sx q[1];
rz(-0.084852858) q[1];
rz(1.3121698) q[3];
sx q[3];
rz(-1.7365518) q[3];
sx q[3];
rz(2.7639619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(-1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(-3.0338874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90751326) q[0];
sx q[0];
rz(-0.36511746) q[0];
sx q[0];
rz(1.3097197) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88476752) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(1.2203072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50385034) q[1];
sx q[1];
rz(-2.3554194) q[1];
sx q[1];
rz(2.1450858) q[1];
rz(2.0884573) q[3];
sx q[3];
rz(-2.5128551) q[3];
sx q[3];
rz(-0.79013463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.084215) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(0.28953141) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(-0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(1.6311197) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(1.0169792) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74422979) q[0];
sx q[0];
rz(-1.1083535) q[0];
sx q[0];
rz(0.37325333) q[0];
rz(-1.594627) q[2];
sx q[2];
rz(-1.7719291) q[2];
sx q[2];
rz(-2.3383274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.14156995) q[1];
sx q[1];
rz(-2.0380028) q[1];
sx q[1];
rz(2.1187374) q[1];
x q[2];
rz(-0.53652699) q[3];
sx q[3];
rz(-1.121472) q[3];
sx q[3];
rz(-2.7262053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54414526) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(-2.5069359) q[2];
rz(2.0641816) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(-0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.2197781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49918136) q[0];
sx q[0];
rz(-1.5326323) q[0];
sx q[0];
rz(0.092060815) q[0];
rz(0.90791038) q[2];
sx q[2];
rz(-2.6320576) q[2];
sx q[2];
rz(-0.31928911) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19983229) q[1];
sx q[1];
rz(-2.4344749) q[1];
sx q[1];
rz(-1.1862399) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7759833) q[3];
sx q[3];
rz(-0.41654166) q[3];
sx q[3];
rz(-0.52683559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62961489) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(-3.1043502) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(1.1720852) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(0.40922871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17657875) q[0];
sx q[0];
rz(-1.8518378) q[0];
sx q[0];
rz(-2.2289508) q[0];
rz(-pi) q[1];
rz(-2.6166603) q[2];
sx q[2];
rz(-1.8443174) q[2];
sx q[2];
rz(-0.33821019) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.65539) q[1];
sx q[1];
rz(-1.4185925) q[1];
sx q[1];
rz(-2.5712719) q[1];
rz(-1.8141361) q[3];
sx q[3];
rz(-1.0662603) q[3];
sx q[3];
rz(-0.53989053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(-3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(-1.508629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2471837) q[0];
sx q[0];
rz(-1.4100473) q[0];
sx q[0];
rz(2.2649691) q[0];
rz(-pi) q[1];
rz(2.1855418) q[2];
sx q[2];
rz(-1.6455257) q[2];
sx q[2];
rz(-1.1751307) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0548045) q[1];
sx q[1];
rz(-2.0874223) q[1];
sx q[1];
rz(0.36887849) q[1];
rz(-2.2739201) q[3];
sx q[3];
rz(-1.5276067) q[3];
sx q[3];
rz(-1.6004576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(0.56525266) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(1.7601097) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0488659) q[0];
sx q[0];
rz(-0.5038358) q[0];
sx q[0];
rz(0.086765246) q[0];
rz(2.6673615) q[2];
sx q[2];
rz(-1.4719163) q[2];
sx q[2];
rz(-0.41254166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5188462) q[1];
sx q[1];
rz(-1.985605) q[1];
sx q[1];
rz(-2.3153789) q[1];
rz(-pi) q[2];
rz(-0.28678203) q[3];
sx q[3];
rz(-2.2638595) q[3];
sx q[3];
rz(-1.1657438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.4505724) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(-0.59664574) q[2];
rz(-2.6560442) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.8687517) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(-1.7455208) q[2];
sx q[2];
rz(-0.25824418) q[2];
sx q[2];
rz(1.2373409) q[2];
rz(-2.6315401) q[3];
sx q[3];
rz(-1.4997713) q[3];
sx q[3];
rz(-0.98239519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
