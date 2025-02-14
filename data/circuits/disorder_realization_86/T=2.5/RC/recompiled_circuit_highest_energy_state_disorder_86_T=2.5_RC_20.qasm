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
rz(0.64492172) q[0];
sx q[0];
rz(2.6725197) q[0];
sx q[0];
rz(7.2735431) q[0];
rz(-1.1733836) q[1];
sx q[1];
rz(-0.85964179) q[1];
sx q[1];
rz(0.8492066) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.828091) q[0];
sx q[0];
rz(-1.573631) q[0];
sx q[0];
rz(1.536973) q[0];
rz(2.8013632) q[2];
sx q[2];
rz(-0.89982596) q[2];
sx q[2];
rz(0.036341993) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4150923) q[1];
sx q[1];
rz(-2.1884568) q[1];
sx q[1];
rz(-2.8128977) q[1];
rz(1.3732713) q[3];
sx q[3];
rz(-0.53086262) q[3];
sx q[3];
rz(1.6429176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90521705) q[2];
sx q[2];
rz(-1.507501) q[2];
sx q[2];
rz(2.6017453) q[2];
rz(-1.5652462) q[3];
sx q[3];
rz(-0.40894517) q[3];
sx q[3];
rz(1.7319771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09963116) q[0];
sx q[0];
rz(-2.4653682) q[0];
sx q[0];
rz(-1.3270295) q[0];
rz(2.3565893) q[1];
sx q[1];
rz(-1.4229341) q[1];
sx q[1];
rz(0.50055093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091565) q[0];
sx q[0];
rz(-1.5958428) q[0];
sx q[0];
rz(-0.5513726) q[0];
rz(1.2039866) q[2];
sx q[2];
rz(-1.7483091) q[2];
sx q[2];
rz(2.8457763) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57823955) q[1];
sx q[1];
rz(-2.1643442) q[1];
sx q[1];
rz(-0.4078354) q[1];
rz(-0.80970069) q[3];
sx q[3];
rz(-1.9040344) q[3];
sx q[3];
rz(-2.7908195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0565722) q[2];
sx q[2];
rz(-0.73740059) q[2];
sx q[2];
rz(-2.6596587) q[2];
rz(-0.36007544) q[3];
sx q[3];
rz(-1.998338) q[3];
sx q[3];
rz(-0.20865194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29927403) q[0];
sx q[0];
rz(-1.1085008) q[0];
sx q[0];
rz(-2.6639248) q[0];
rz(2.281669) q[1];
sx q[1];
rz(-1.0483024) q[1];
sx q[1];
rz(-2.8544676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5764113) q[0];
sx q[0];
rz(-0.14248304) q[0];
sx q[0];
rz(1.0866685) q[0];
x q[1];
rz(2.6052104) q[2];
sx q[2];
rz(-2.7349758) q[2];
sx q[2];
rz(2.0970203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7487919) q[1];
sx q[1];
rz(-2.254389) q[1];
sx q[1];
rz(-3.1303065) q[1];
rz(-pi) q[2];
rz(1.1895387) q[3];
sx q[3];
rz(-1.4917456) q[3];
sx q[3];
rz(1.4850657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75256789) q[2];
sx q[2];
rz(-1.3501046) q[2];
sx q[2];
rz(3.1217421) q[2];
rz(-0.94414532) q[3];
sx q[3];
rz(-0.70725924) q[3];
sx q[3];
rz(2.7437362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860564) q[0];
sx q[0];
rz(-2.5113386) q[0];
sx q[0];
rz(1.8408884) q[0];
rz(2.6553254) q[1];
sx q[1];
rz(-1.9673037) q[1];
sx q[1];
rz(0.035331443) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080950532) q[0];
sx q[0];
rz(-1.544702) q[0];
sx q[0];
rz(-1.1423201) q[0];
rz(-pi) q[1];
rz(-0.51574083) q[2];
sx q[2];
rz(-2.6922142) q[2];
sx q[2];
rz(0.55381227) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7481193) q[1];
sx q[1];
rz(-1.5941718) q[1];
sx q[1];
rz(-1.5651902) q[1];
rz(-pi) q[2];
rz(-2.1256623) q[3];
sx q[3];
rz(-2.1240222) q[3];
sx q[3];
rz(-0.31496668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47762927) q[2];
sx q[2];
rz(-0.56325459) q[2];
sx q[2];
rz(-2.8832054) q[2];
rz(0.7102617) q[3];
sx q[3];
rz(-2.5370772) q[3];
sx q[3];
rz(-1.5822423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6579987) q[0];
sx q[0];
rz(-0.76376629) q[0];
sx q[0];
rz(-2.4972231) q[0];
rz(2.1517892) q[1];
sx q[1];
rz(-2.3376696) q[1];
sx q[1];
rz(-2.03233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.155494) q[0];
sx q[0];
rz(-1.3751831) q[0];
sx q[0];
rz(1.0367111) q[0];
rz(-pi) q[1];
rz(2.0844056) q[2];
sx q[2];
rz(-1.6053146) q[2];
sx q[2];
rz(0.90189722) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7444744) q[1];
sx q[1];
rz(-2.3496685) q[1];
sx q[1];
rz(2.9379815) q[1];
rz(-pi) q[2];
rz(2.3886834) q[3];
sx q[3];
rz(-1.782182) q[3];
sx q[3];
rz(1.1803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7295419) q[2];
sx q[2];
rz(-1.5645626) q[2];
sx q[2];
rz(-0.61005074) q[2];
rz(0.080502056) q[3];
sx q[3];
rz(-0.1463612) q[3];
sx q[3];
rz(0.0065053594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7215111) q[0];
sx q[0];
rz(-0.93669909) q[0];
sx q[0];
rz(1.4790081) q[0];
rz(-1.1889907) q[1];
sx q[1];
rz(-0.34591302) q[1];
sx q[1];
rz(1.6993274) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8540031) q[0];
sx q[0];
rz(-0.95335273) q[0];
sx q[0];
rz(2.94728) q[0];
rz(-pi) q[1];
rz(1.5046607) q[2];
sx q[2];
rz(-2.1512262) q[2];
sx q[2];
rz(2.5709267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.14076041) q[1];
sx q[1];
rz(-0.78703431) q[1];
sx q[1];
rz(-2.0443022) q[1];
rz(-pi) q[2];
x q[2];
rz(2.598212) q[3];
sx q[3];
rz(-0.60688775) q[3];
sx q[3];
rz(0.010963765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46001616) q[2];
sx q[2];
rz(-0.67395335) q[2];
sx q[2];
rz(-0.67503929) q[2];
rz(-2.8071844) q[3];
sx q[3];
rz(-1.6655917) q[3];
sx q[3];
rz(-2.7736751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6595031) q[0];
sx q[0];
rz(-2.1559494) q[0];
sx q[0];
rz(3.0159045) q[0];
rz(-0.19371678) q[1];
sx q[1];
rz(-2.2592762) q[1];
sx q[1];
rz(-1.9901336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4070258) q[0];
sx q[0];
rz(-2.1687963) q[0];
sx q[0];
rz(2.1728188) q[0];
rz(-pi) q[1];
rz(2.861768) q[2];
sx q[2];
rz(-1.9283659) q[2];
sx q[2];
rz(-1.8574865) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0380504) q[1];
sx q[1];
rz(-1.0443496) q[1];
sx q[1];
rz(1.0564694) q[1];
x q[2];
rz(-2.4667593) q[3];
sx q[3];
rz(-1.1911567) q[3];
sx q[3];
rz(1.4745431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2842399) q[2];
sx q[2];
rz(-1.7809296) q[2];
sx q[2];
rz(0.24687684) q[2];
rz(-2.2829368) q[3];
sx q[3];
rz(-2.1571428) q[3];
sx q[3];
rz(0.27913678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7849279) q[0];
sx q[0];
rz(-2.6981638) q[0];
sx q[0];
rz(0.32522935) q[0];
rz(-0.52109703) q[1];
sx q[1];
rz(-2.5083713) q[1];
sx q[1];
rz(0.31347832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26599182) q[0];
sx q[0];
rz(-0.82339215) q[0];
sx q[0];
rz(0.70988795) q[0];
rz(-pi) q[1];
rz(1.7960848) q[2];
sx q[2];
rz(-1.6367925) q[2];
sx q[2];
rz(0.57541945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65113803) q[1];
sx q[1];
rz(-0.50781194) q[1];
sx q[1];
rz(2.3976624) q[1];
rz(-2.0920803) q[3];
sx q[3];
rz(-1.7465599) q[3];
sx q[3];
rz(2.5960032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7079805) q[2];
sx q[2];
rz(-1.2182451) q[2];
sx q[2];
rz(-1.2797959) q[2];
rz(-2.4158939) q[3];
sx q[3];
rz(-1.9819219) q[3];
sx q[3];
rz(-0.80997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5028266) q[0];
sx q[0];
rz(-1.1932729) q[0];
sx q[0];
rz(0.1499114) q[0];
rz(-1.2431078) q[1];
sx q[1];
rz(-1.5877692) q[1];
sx q[1];
rz(-1.4814203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4325162) q[0];
sx q[0];
rz(-3.0750599) q[0];
sx q[0];
rz(-2.3155022) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3636087) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(-0.65226511) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2592377) q[1];
sx q[1];
rz(-1.9505525) q[1];
sx q[1];
rz(1.3391375) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8152913) q[3];
sx q[3];
rz(-1.4473379) q[3];
sx q[3];
rz(-2.3924912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3159065) q[2];
sx q[2];
rz(-1.76182) q[2];
sx q[2];
rz(-2.1659577) q[2];
rz(-1.1286831) q[3];
sx q[3];
rz(-2.5895139) q[3];
sx q[3];
rz(0.25477195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16633701) q[0];
sx q[0];
rz(-1.0003426) q[0];
sx q[0];
rz(0.16217232) q[0];
rz(-1.8099248) q[1];
sx q[1];
rz(-0.63260308) q[1];
sx q[1];
rz(0.046028927) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9105658) q[0];
sx q[0];
rz(-1.5648769) q[0];
sx q[0];
rz(1.5886515) q[0];
rz(2.677605) q[2];
sx q[2];
rz(-1.2210238) q[2];
sx q[2];
rz(-2.2217563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0303453) q[1];
sx q[1];
rz(-1.1702036) q[1];
sx q[1];
rz(0.16040332) q[1];
rz(-1.7415984) q[3];
sx q[3];
rz(-1.2360337) q[3];
sx q[3];
rz(2.6790909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33596805) q[2];
sx q[2];
rz(-0.38463548) q[2];
sx q[2];
rz(1.1495205) q[2];
rz(2.6286821) q[3];
sx q[3];
rz(-1.7549113) q[3];
sx q[3];
rz(2.8368867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.649986) q[0];
sx q[0];
rz(-0.91616022) q[0];
sx q[0];
rz(1.1471163) q[0];
rz(0.06123771) q[1];
sx q[1];
rz(-1.9495268) q[1];
sx q[1];
rz(-1.3757642) q[1];
rz(0.55228615) q[2];
sx q[2];
rz(-0.88197642) q[2];
sx q[2];
rz(-2.3563202) q[2];
rz(-1.9763532) q[3];
sx q[3];
rz(-0.83409393) q[3];
sx q[3];
rz(0.50611151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
