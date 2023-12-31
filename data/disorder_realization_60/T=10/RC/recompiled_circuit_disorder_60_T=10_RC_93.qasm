OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1922167) q[0];
sx q[0];
rz(-2.0944216) q[0];
sx q[0];
rz(-0.068724364) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(1.2083763) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87586227) q[0];
sx q[0];
rz(-3.0348572) q[0];
sx q[0];
rz(-2.021832) q[0];
rz(-3/(10*pi)) q[2];
sx q[2];
rz(-3.0152233) q[2];
sx q[2];
rz(-2.7356844) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2627416) q[1];
sx q[1];
rz(-2.1858715) q[1];
sx q[1];
rz(-1.7064852) q[1];
rz(-2.4083706) q[3];
sx q[3];
rz(-0.65922046) q[3];
sx q[3];
rz(0.99430195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8157114) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(0.74716032) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0242457) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(2.970447) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(-0.29719621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3127808) q[0];
sx q[0];
rz(-3.045407) q[0];
sx q[0];
rz(1.0440774) q[0];
rz(-0.53571312) q[2];
sx q[2];
rz(-0.78768724) q[2];
sx q[2];
rz(0.96131575) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52642194) q[1];
sx q[1];
rz(-0.63262683) q[1];
sx q[1];
rz(0.97140615) q[1];
x q[2];
rz(0.76241775) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(3.1114515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16195665) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(-0.61398181) q[2];
rz(-0.87614122) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(0.63703018) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-2.8339548) q[0];
rz(0.74854198) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(-0.83980733) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6162286) q[0];
sx q[0];
rz(-1.9814241) q[0];
sx q[0];
rz(1.8812268) q[0];
rz(-pi) q[1];
rz(0.33939056) q[2];
sx q[2];
rz(-1.6600142) q[2];
sx q[2];
rz(0.81775507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0864598) q[1];
sx q[1];
rz(-1.8756525) q[1];
sx q[1];
rz(-2.6677368) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8057459) q[3];
sx q[3];
rz(-0.64628212) q[3];
sx q[3];
rz(2.7502053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.86747375) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(-1.9258202) q[3];
sx q[3];
rz(-0.35651818) q[3];
sx q[3];
rz(1.5095476) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76386219) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(2.6960301) q[0];
rz(2.5207649) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-0.96558085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8646116) q[0];
sx q[0];
rz(-0.81482139) q[0];
sx q[0];
rz(-2.7576202) q[0];
rz(1.5358652) q[2];
sx q[2];
rz(-2.4128649) q[2];
sx q[2];
rz(-0.59128371) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4703776) q[1];
sx q[1];
rz(-2.885985) q[1];
sx q[1];
rz(1.4443936) q[1];
x q[2];
rz(0.86430092) q[3];
sx q[3];
rz(-0.14926499) q[3];
sx q[3];
rz(0.60993689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-0.92932534) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716361) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(1.2840282) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(-1.0481542) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3583583) q[0];
sx q[0];
rz(-2.2336707) q[0];
sx q[0];
rz(0.62028424) q[0];
rz(-2.3088147) q[2];
sx q[2];
rz(-1.1156429) q[2];
sx q[2];
rz(-0.27311329) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.594164) q[1];
sx q[1];
rz(-2.8058726) q[1];
sx q[1];
rz(2.8572542) q[1];
x q[2];
rz(-0.9976451) q[3];
sx q[3];
rz(-2.3754658) q[3];
sx q[3];
rz(2.0875974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(-2.0488996) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(0.95935217) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(-0.59610468) q[0];
rz(1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.2449107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0743474) q[0];
sx q[0];
rz(-2.0739569) q[0];
sx q[0];
rz(0.36672451) q[0];
rz(-pi) q[1];
rz(-2.6558098) q[2];
sx q[2];
rz(-1.2885639) q[2];
sx q[2];
rz(0.9529875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5925496) q[1];
sx q[1];
rz(-2.2144496) q[1];
sx q[1];
rz(0.93306577) q[1];
x q[2];
rz(1.3390433) q[3];
sx q[3];
rz(-1.8125121) q[3];
sx q[3];
rz(1.7237323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(-0.088430017) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769619) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(-2.8616469) q[0];
rz(-1.6784558) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(0.25269145) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1278909) q[0];
sx q[0];
rz(-1.863592) q[0];
sx q[0];
rz(-3.091759) q[0];
x q[1];
rz(-1.5964609) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(0.54021013) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.686324) q[1];
sx q[1];
rz(-0.9372006) q[1];
sx q[1];
rz(2.2151161) q[1];
x q[2];
rz(2.2884376) q[3];
sx q[3];
rz(-2.1907638) q[3];
sx q[3];
rz(-0.48228797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(-1.6097216) q[2];
rz(-1.948471) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(-2.8727818) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(2.862646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69142197) q[0];
sx q[0];
rz(-2.2922463) q[0];
sx q[0];
rz(0.88766092) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0729162) q[2];
sx q[2];
rz(-1.4197822) q[2];
sx q[2];
rz(-0.47889027) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8262417) q[1];
sx q[1];
rz(-2.3126174) q[1];
sx q[1];
rz(-2.6226603) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14848498) q[3];
sx q[3];
rz(-1.0035702) q[3];
sx q[3];
rz(2.4384769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3387317) q[2];
sx q[2];
rz(-1.4638476) q[2];
sx q[2];
rz(-1.5489244) q[2];
rz(-2.0907949) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(1.2532225) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(-1.1463096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.779594) q[0];
sx q[0];
rz(-0.75095526) q[0];
sx q[0];
rz(0.10303084) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65593221) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(2.7149372) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6625001) q[1];
sx q[1];
rz(-0.88516419) q[1];
sx q[1];
rz(-0.094866026) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54543145) q[3];
sx q[3];
rz(-2.430393) q[3];
sx q[3];
rz(-1.7033073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(-2.1257607) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(-2.773496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(3.1273499) q[0];
rz(0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(1.7600118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.325557) q[0];
sx q[0];
rz(-2.9986458) q[0];
sx q[0];
rz(-1.4965034) q[0];
rz(-pi) q[1];
rz(3.0773452) q[2];
sx q[2];
rz(-2.148743) q[2];
sx q[2];
rz(2.3648175) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.7144202) q[1];
sx q[1];
rz(2.1302057) q[1];
rz(0.53346177) q[3];
sx q[3];
rz(-1.7603612) q[3];
sx q[3];
rz(1.5272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(-2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(-0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0902696) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(-1.8016215) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(-2.6132265) q[2];
sx q[2];
rz(-0.7893723) q[2];
sx q[2];
rz(-1.5783725) q[2];
rz(2.1918478) q[3];
sx q[3];
rz(-1.2093778) q[3];
sx q[3];
rz(-2.9686684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
