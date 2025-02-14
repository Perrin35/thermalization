OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1005062) q[0];
sx q[0];
rz(-2.9334928) q[0];
sx q[0];
rz(-2.7383374) q[0];
rz(-0.23303214) q[1];
sx q[1];
rz(-1.4401399) q[1];
sx q[1];
rz(-2.9177102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1772645) q[0];
sx q[0];
rz(-1.7814264) q[0];
sx q[0];
rz(0.62646477) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12822784) q[2];
sx q[2];
rz(-2.2710481) q[2];
sx q[2];
rz(-2.4403404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2157872) q[1];
sx q[1];
rz(-0.69725446) q[1];
sx q[1];
rz(1.2165192) q[1];
x q[2];
rz(2.3305064) q[3];
sx q[3];
rz(-1.038365) q[3];
sx q[3];
rz(0.48871751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69748059) q[2];
sx q[2];
rz(-0.31284249) q[2];
sx q[2];
rz(0.035813896) q[2];
rz(-0.73016417) q[3];
sx q[3];
rz(-1.1532447) q[3];
sx q[3];
rz(-2.9353976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.25074729) q[0];
sx q[0];
rz(-2.146281) q[0];
sx q[0];
rz(0.19749755) q[0];
rz(-0.47710553) q[1];
sx q[1];
rz(-0.66480607) q[1];
sx q[1];
rz(-0.8055996) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6960736) q[0];
sx q[0];
rz(-1.7564) q[0];
sx q[0];
rz(0.66849065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7609772) q[2];
sx q[2];
rz(-1.5089436) q[2];
sx q[2];
rz(-2.672489) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9220613) q[1];
sx q[1];
rz(-2.1967255) q[1];
sx q[1];
rz(-2.3971167) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74975852) q[3];
sx q[3];
rz(-2.9364412) q[3];
sx q[3];
rz(-1.6501901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8311367) q[2];
sx q[2];
rz(-1.6418991) q[2];
sx q[2];
rz(-1.4031225) q[2];
rz(2.5095615) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(3.0145751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8359351) q[0];
sx q[0];
rz(-0.62863612) q[0];
sx q[0];
rz(1.0868616) q[0];
rz(-2.944259) q[1];
sx q[1];
rz(-1.6355762) q[1];
sx q[1];
rz(-0.33263439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5088288) q[0];
sx q[0];
rz(-0.45464215) q[0];
sx q[0];
rz(-1.8931382) q[0];
rz(-3.0789571) q[2];
sx q[2];
rz(-0.83683678) q[2];
sx q[2];
rz(-2.1750185) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2283295) q[1];
sx q[1];
rz(-2.8614534) q[1];
sx q[1];
rz(0.88714169) q[1];
rz(-pi) q[2];
rz(-2.8841697) q[3];
sx q[3];
rz(-0.74353131) q[3];
sx q[3];
rz(-0.0071255077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3745594) q[2];
sx q[2];
rz(-1.5908396) q[2];
sx q[2];
rz(1.4683051) q[2];
rz(3.1324006) q[3];
sx q[3];
rz(-2.6720948) q[3];
sx q[3];
rz(-0.79206842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.511635) q[0];
sx q[0];
rz(-2.503643) q[0];
sx q[0];
rz(1.4427503) q[0];
rz(1.0244145) q[1];
sx q[1];
rz(-1.9099648) q[1];
sx q[1];
rz(-2.3146497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8762739) q[0];
sx q[0];
rz(-1.4965222) q[0];
sx q[0];
rz(1.386388) q[0];
x q[1];
rz(-0.154279) q[2];
sx q[2];
rz(-0.88866975) q[2];
sx q[2];
rz(0.23274225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1093028) q[1];
sx q[1];
rz(-2.656611) q[1];
sx q[1];
rz(-0.46320199) q[1];
rz(-0.8884646) q[3];
sx q[3];
rz(-1.6649705) q[3];
sx q[3];
rz(1.7674131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8168762) q[2];
sx q[2];
rz(-2.3574895) q[2];
sx q[2];
rz(-1.2775705) q[2];
rz(2.4534524) q[3];
sx q[3];
rz(-1.0253996) q[3];
sx q[3];
rz(-0.96243206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14030305) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(-2.9241614) q[0];
rz(-0.28542074) q[1];
sx q[1];
rz(-2.5358584) q[1];
sx q[1];
rz(0.4981471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0131684) q[0];
sx q[0];
rz(-1.2449055) q[0];
sx q[0];
rz(-2.2504435) q[0];
rz(-pi) q[1];
rz(-0.26136036) q[2];
sx q[2];
rz(-2.0730505) q[2];
sx q[2];
rz(2.4136191) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0302311) q[1];
sx q[1];
rz(-1.7788789) q[1];
sx q[1];
rz(1.0743121) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6400233) q[3];
sx q[3];
rz(-0.49276982) q[3];
sx q[3];
rz(1.2214965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.51776) q[2];
sx q[2];
rz(-1.0972923) q[2];
sx q[2];
rz(0.1499873) q[2];
rz(-0.013966694) q[3];
sx q[3];
rz(-1.5499127) q[3];
sx q[3];
rz(-0.51378957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496395) q[0];
sx q[0];
rz(-0.54170251) q[0];
sx q[0];
rz(-2.3288222) q[0];
rz(-1.4179519) q[1];
sx q[1];
rz(-1.6287454) q[1];
sx q[1];
rz(-2.0779804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2246216) q[0];
sx q[0];
rz(-1.8993036) q[0];
sx q[0];
rz(-0.55060951) q[0];
rz(-pi) q[1];
rz(3.0399357) q[2];
sx q[2];
rz(-0.16451193) q[2];
sx q[2];
rz(0.042334231) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2460766) q[1];
sx q[1];
rz(-2.5427596) q[1];
sx q[1];
rz(-0.020222874) q[1];
x q[2];
rz(1.316458) q[3];
sx q[3];
rz(-1.3876378) q[3];
sx q[3];
rz(1.1381799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26760462) q[2];
sx q[2];
rz(-0.559811) q[2];
sx q[2];
rz(-1.416729) q[2];
rz(-0.060976107) q[3];
sx q[3];
rz(-0.91106001) q[3];
sx q[3];
rz(1.4643668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7199719) q[0];
sx q[0];
rz(-2.3201729) q[0];
sx q[0];
rz(0.11909568) q[0];
rz(-0.47856092) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(2.4518769) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091662377) q[0];
sx q[0];
rz(-2.1436467) q[0];
sx q[0];
rz(-0.86688231) q[0];
rz(-pi) q[1];
rz(1.941626) q[2];
sx q[2];
rz(-2.0304108) q[2];
sx q[2];
rz(0.72061611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7642245) q[1];
sx q[1];
rz(-0.86478327) q[1];
sx q[1];
rz(3.0824667) q[1];
rz(-pi) q[2];
rz(2.2834217) q[3];
sx q[3];
rz(-1.6827876) q[3];
sx q[3];
rz(-1.5369161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0927642) q[2];
sx q[2];
rz(-0.060828716) q[2];
sx q[2];
rz(2.0737341) q[2];
rz(2.974406) q[3];
sx q[3];
rz(-1.7696295) q[3];
sx q[3];
rz(-0.13944496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6086513) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(-0.90840489) q[0];
rz(2.1482229) q[1];
sx q[1];
rz(-2.9790331) q[1];
sx q[1];
rz(2.2672674) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217978) q[0];
sx q[0];
rz(-0.13331977) q[0];
sx q[0];
rz(-0.13102417) q[0];
rz(-1.3130593) q[2];
sx q[2];
rz(-1.589587) q[2];
sx q[2];
rz(-0.37237871) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5914375) q[1];
sx q[1];
rz(-0.14796013) q[1];
sx q[1];
rz(1.1569389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0197952) q[3];
sx q[3];
rz(-1.8809109) q[3];
sx q[3];
rz(1.9581025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61451644) q[2];
sx q[2];
rz(-0.15915844) q[2];
sx q[2];
rz(0.85421526) q[2];
rz(-2.6863875) q[3];
sx q[3];
rz(-0.56536094) q[3];
sx q[3];
rz(-2.8860886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47827569) q[0];
sx q[0];
rz(-1.217696) q[0];
sx q[0];
rz(1.4772557) q[0];
rz(1.5746337) q[1];
sx q[1];
rz(-2.4576371) q[1];
sx q[1];
rz(-2.4050567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57845302) q[0];
sx q[0];
rz(-1.805465) q[0];
sx q[0];
rz(-2.0521972) q[0];
x q[1];
rz(-2.8868448) q[2];
sx q[2];
rz(-2.7494135) q[2];
sx q[2];
rz(0.90975159) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29727259) q[1];
sx q[1];
rz(-1.2579009) q[1];
sx q[1];
rz(-1.827153) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1677954) q[3];
sx q[3];
rz(-1.6000119) q[3];
sx q[3];
rz(0.67724281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3972724) q[2];
sx q[2];
rz(-0.87817764) q[2];
sx q[2];
rz(1.6142023) q[2];
rz(-2.7496036) q[3];
sx q[3];
rz(-1.9422928) q[3];
sx q[3];
rz(-0.049662445) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.117332) q[0];
sx q[0];
rz(-1.4735104) q[0];
sx q[0];
rz(2.9397553) q[0];
rz(-0.2233389) q[1];
sx q[1];
rz(-0.52450648) q[1];
sx q[1];
rz(2.7490659) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46014547) q[0];
sx q[0];
rz(-2.0325615) q[0];
sx q[0];
rz(0.11790922) q[0];
x q[1];
rz(0.052458737) q[2];
sx q[2];
rz(-2.3067368) q[2];
sx q[2];
rz(0.8142161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10822693) q[1];
sx q[1];
rz(-1.1443597) q[1];
sx q[1];
rz(-0.99797499) q[1];
rz(-0.5397615) q[3];
sx q[3];
rz(-2.5934412) q[3];
sx q[3];
rz(-2.2659311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7212123) q[2];
sx q[2];
rz(-2.1537697) q[2];
sx q[2];
rz(0.58050275) q[2];
rz(-0.37603363) q[3];
sx q[3];
rz(-0.97085634) q[3];
sx q[3];
rz(1.5413126) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15039438) q[0];
sx q[0];
rz(-1.5491485) q[0];
sx q[0];
rz(-1.3843672) q[0];
rz(-1.4906384) q[1];
sx q[1];
rz(-1.8883659) q[1];
sx q[1];
rz(-2.2813003) q[1];
rz(0.91928405) q[2];
sx q[2];
rz(-2.2625661) q[2];
sx q[2];
rz(1.7023466) q[2];
rz(-2.8297791) q[3];
sx q[3];
rz(-1.5494294) q[3];
sx q[3];
rz(2.6904306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
