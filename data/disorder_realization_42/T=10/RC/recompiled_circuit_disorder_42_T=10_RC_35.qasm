OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(-2.7789814) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(-0.4904823) q[1];
sx q[1];
rz(2.7999556) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65389079) q[0];
sx q[0];
rz(-2.747295) q[0];
sx q[0];
rz(2.8173692) q[0];
x q[1];
rz(-0.59092893) q[2];
sx q[2];
rz(-1.8701631) q[2];
sx q[2];
rz(-0.17609827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7340378) q[1];
sx q[1];
rz(-1.8681548) q[1];
sx q[1];
rz(2.107891) q[1];
rz(-pi) q[2];
rz(-1.0757252) q[3];
sx q[3];
rz(-0.44701156) q[3];
sx q[3];
rz(1.8398374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59387702) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(0.95648742) q[2];
rz(-2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.1454426) q[0];
rz(-1.068813) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(2.4157445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024591) q[0];
sx q[0];
rz(-2.4665138) q[0];
sx q[0];
rz(-2.549987) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1109353) q[2];
sx q[2];
rz(-1.0265372) q[2];
sx q[2];
rz(-2.0366675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5686381) q[1];
sx q[1];
rz(-1.2619839) q[1];
sx q[1];
rz(3.1152578) q[1];
rz(-1.3594567) q[3];
sx q[3];
rz(-0.95892116) q[3];
sx q[3];
rz(-2.8733727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-0.99728161) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230351) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(1.1285271) q[0];
rz(1.3035125) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(0.9054786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826686) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(-0.24567901) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0037903) q[2];
sx q[2];
rz(-0.86647063) q[2];
sx q[2];
rz(1.1906884) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8454051) q[1];
sx q[1];
rz(-2.2803109) q[1];
sx q[1];
rz(0.097464949) q[1];
rz(-2.2842992) q[3];
sx q[3];
rz(-1.6033353) q[3];
sx q[3];
rz(-1.9705704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.7604312) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(0.9872438) q[2];
rz(0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-0.12606829) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(-1.1674081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084598736) q[0];
sx q[0];
rz(-1.7718678) q[0];
sx q[0];
rz(2.3720471) q[0];
rz(-1.7502968) q[2];
sx q[2];
rz(-1.0738392) q[2];
sx q[2];
rz(-2.3617982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3525225) q[1];
sx q[1];
rz(-1.8348502) q[1];
sx q[1];
rz(-2.4697484) q[1];
rz(-pi) q[2];
rz(-0.64658029) q[3];
sx q[3];
rz(-2.1255891) q[3];
sx q[3];
rz(-0.085689714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(0.63151044) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773848) q[0];
sx q[0];
rz(-3.003484) q[0];
sx q[0];
rz(2.6389129) q[0];
rz(2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(0.46708333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57468092) q[0];
sx q[0];
rz(-1.7705288) q[0];
sx q[0];
rz(1.9499705) q[0];
rz(0.58606834) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(0.52473247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5284164) q[1];
sx q[1];
rz(-1.3532552) q[1];
sx q[1];
rz(1.1385285) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2561982) q[3];
sx q[3];
rz(-1.6939031) q[3];
sx q[3];
rz(0.76243329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6749394) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(-2.6494027) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(0.870134) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(-1.4280691) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8498358) q[0];
sx q[0];
rz(-1.0514326) q[0];
sx q[0];
rz(0.22626466) q[0];
x q[1];
rz(1.0322083) q[2];
sx q[2];
rz(-1.7787873) q[2];
sx q[2];
rz(-1.6230735) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9073346) q[1];
sx q[1];
rz(-1.9252216) q[1];
sx q[1];
rz(-0.025269421) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1268678) q[3];
sx q[3];
rz(-2.428599) q[3];
sx q[3];
rz(1.4006437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(0.39624828) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8783766) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(0.011627442) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41650018) q[0];
sx q[0];
rz(-1.9511127) q[0];
sx q[0];
rz(-1.8561383) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3064012) q[2];
sx q[2];
rz(-3.1036658) q[2];
sx q[2];
rz(1.7131625) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7203622) q[1];
sx q[1];
rz(-2.1167397) q[1];
sx q[1];
rz(0.027688428) q[1];
rz(2.4680448) q[3];
sx q[3];
rz(-1.3140956) q[3];
sx q[3];
rz(1.4207157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(-0.28798506) q[2];
rz(1.9912432) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(-1.9294552) q[0];
rz(2.7638226) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(1.6479962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15570116) q[0];
sx q[0];
rz(-2.771286) q[0];
sx q[0];
rz(0.62472384) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0817501) q[2];
sx q[2];
rz(-1.7027731) q[2];
sx q[2];
rz(-3.0795385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89477506) q[1];
sx q[1];
rz(-0.90118876) q[1];
sx q[1];
rz(-2.0991904) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11843189) q[3];
sx q[3];
rz(-1.0716972) q[3];
sx q[3];
rz(-1.1008319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3528072) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(-1.8898233) q[2];
rz(2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-0.10929508) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733474) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(0.1782724) q[0];
rz(-0.072862236) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(-1.1791621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68547738) q[0];
sx q[0];
rz(-0.59887409) q[0];
sx q[0];
rz(2.6045538) q[0];
rz(-0.91782848) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(2.3248621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0852016) q[1];
sx q[1];
rz(-1.3136275) q[1];
sx q[1];
rz(-2.609054) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0719028) q[3];
sx q[3];
rz(-1.9074829) q[3];
sx q[3];
rz(-1.5130373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4420085) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(2.1098095) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(2.9163196) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28829065) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(0.061696079) q[0];
rz(1.306698) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(2.0526989) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63554791) q[0];
sx q[0];
rz(-1.5061) q[0];
sx q[0];
rz(-0.062775469) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(2.0053787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.006146487) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(-2.1250493) q[1];
rz(-2.0652566) q[3];
sx q[3];
rz(-2.7571602) q[3];
sx q[3];
rz(-1.2151375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5237727) q[2];
sx q[2];
rz(-0.68796316) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(-0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4927647) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(-1.6246673) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(0.4060612) q[2];
sx q[2];
rz(-2.1880697) q[2];
sx q[2];
rz(-0.38068117) q[2];
rz(1.6886961) q[3];
sx q[3];
rz(-2.8040734) q[3];
sx q[3];
rz(0.87344195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
