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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87586227) q[0];
sx q[0];
rz(-3.0348572) q[0];
sx q[0];
rz(1.1197607) q[0];
rz(-0.12579972) q[2];
sx q[2];
rz(-1.5828136) q[2];
sx q[2];
rz(-1.2596241) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.61332834) q[1];
sx q[1];
rz(-1.6815038) q[1];
sx q[1];
rz(0.61943357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52238676) q[3];
sx q[3];
rz(-1.99317) q[3];
sx q[3];
rz(-1.9463604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.6750083) q[2];
rz(-2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(-0.74716032) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0242457) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(1.9673989) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(-0.29719621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3127808) q[0];
sx q[0];
rz(-0.096185616) q[0];
sx q[0];
rz(1.0440774) q[0];
rz(-pi) q[1];
rz(-2.6058795) q[2];
sx q[2];
rz(-2.3539054) q[2];
sx q[2];
rz(-2.1802769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5479996) q[1];
sx q[1];
rz(-1.2307234) q[1];
sx q[1];
rz(1.026457) q[1];
x q[2];
rz(-0.76241775) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(-0.61398181) q[2];
rz(-2.2654514) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0456332) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(2.8339548) q[0];
rz(2.3930507) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(2.3017853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6162286) q[0];
sx q[0];
rz(-1.1601686) q[0];
sx q[0];
rz(-1.2603659) q[0];
rz(-pi) q[1];
rz(-0.26250458) q[2];
sx q[2];
rz(-2.7911107) q[2];
sx q[2];
rz(-2.1413435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.055132853) q[1];
sx q[1];
rz(-1.8756525) q[1];
sx q[1];
rz(2.6677368) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3358467) q[3];
sx q[3];
rz(-0.64628212) q[3];
sx q[3];
rz(0.39138734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(-1.3179368) q[2];
rz(-1.2157724) q[3];
sx q[3];
rz(-0.35651818) q[3];
sx q[3];
rz(1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76386219) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(-0.44556251) q[0];
rz(0.62082779) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(2.1760118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5774028) q[0];
sx q[0];
rz(-1.2947384) q[0];
sx q[0];
rz(-2.3645556) q[0];
x q[1];
rz(1.6057274) q[2];
sx q[2];
rz(-2.4128649) q[2];
sx q[2];
rz(-2.5503089) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.671215) q[1];
sx q[1];
rz(-0.2556076) q[1];
sx q[1];
rz(-1.4443936) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2772917) q[3];
sx q[3];
rz(-2.9923277) q[3];
sx q[3];
rz(-2.5316558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(0.6435414) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(2.1381366) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716361) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(-1.8575645) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(2.0934385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6266039) q[0];
sx q[0];
rz(-2.0467313) q[0];
sx q[0];
rz(-2.3355052) q[0];
rz(0.58446144) q[2];
sx q[2];
rz(-0.9218773) q[2];
sx q[2];
rz(1.6780168) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2939799) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(-1.4732248) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6606779) q[3];
sx q[3];
rz(-0.9489343) q[3];
sx q[3];
rz(2.8180168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2709048) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(2.2407545) q[2];
rz(-2.0488996) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(-1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95935217) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(-0.59610468) q[0];
rz(-1.4959363) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.8966819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67958528) q[0];
sx q[0];
rz(-1.2512659) q[0];
sx q[0];
rz(-1.0380448) q[0];
rz(-pi) q[1];
rz(-1.8876569) q[2];
sx q[2];
rz(-2.0358026) q[2];
sx q[2];
rz(0.76380619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.54904304) q[1];
sx q[1];
rz(-0.92714308) q[1];
sx q[1];
rz(0.93306577) q[1];
x q[2];
rz(-0.75002589) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(-0.94543524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.231455) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(0.088430017) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(-2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6769619) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(-2.8616469) q[0];
rz(-1.4631368) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(2.8889012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1278909) q[0];
sx q[0];
rz(-1.863592) q[0];
sx q[0];
rz(-0.049833628) q[0];
rz(-1.5451317) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(2.6013825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44811571) q[1];
sx q[1];
rz(-0.87065334) q[1];
sx q[1];
rz(-2.4561988) q[1];
rz(-0.74440646) q[3];
sx q[3];
rz(-0.91074569) q[3];
sx q[3];
rz(-0.50124121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(-1.5318711) q[2];
rz(-1.1931217) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10483345) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(-1.4861134) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(-0.2789467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69142197) q[0];
sx q[0];
rz(-2.2922463) q[0];
sx q[0];
rz(2.2539317) q[0];
rz(1.8770304) q[2];
sx q[2];
rz(-0.5224723) q[2];
sx q[2];
rz(-0.82447169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3153509) q[1];
sx q[1];
rz(-2.3126174) q[1];
sx q[1];
rz(-2.6226603) q[1];
x q[2];
rz(2.9931077) q[3];
sx q[3];
rz(-2.1380224) q[3];
sx q[3];
rz(2.4384769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(1.1999493) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46362296) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(-1.8883702) q[0];
rz(-2.5190917) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.995283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86664591) q[0];
sx q[0];
rz(-1.5005611) q[0];
sx q[0];
rz(0.74830351) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4856604) q[2];
sx q[2];
rz(-2.1428875) q[2];
sx q[2];
rz(2.7149372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1100626) q[1];
sx q[1];
rz(-1.6441802) q[1];
sx q[1];
rz(-2.258638) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5961612) q[3];
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
rz(-1.035707) q[2];
sx q[2];
rz(-1.0158319) q[2];
rz(-0.9097957) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0125473) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(-0.01424271) q[0];
rz(2.2968538) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(1.7600118) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81603564) q[0];
sx q[0];
rz(-0.14294681) q[0];
sx q[0];
rz(-1.6450892) q[0];
rz(1.6689156) q[2];
sx q[2];
rz(-0.58110229) q[2];
sx q[2];
rz(0.89400089) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8091619) q[1];
sx q[1];
rz(-0.5756439) q[1];
sx q[1];
rz(-1.3047421) q[1];
x q[2];
rz(2.6081309) q[3];
sx q[3];
rz(-1.7603612) q[3];
sx q[3];
rz(-1.5272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0766729) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(-2.1515965) q[2];
rz(-0.89407095) q[3];
sx q[3];
rz(-0.48864135) q[3];
sx q[3];
rz(2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0513231) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(-1.3399711) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(-0.71628911) q[2];
sx q[2];
rz(-1.9367957) q[2];
sx q[2];
rz(-0.39762485) q[2];
rz(0.94974489) q[3];
sx q[3];
rz(-1.9322149) q[3];
sx q[3];
rz(0.17292427) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
