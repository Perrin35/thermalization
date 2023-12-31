OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.94937593) q[0];
sx q[0];
rz(-1.047171) q[0];
sx q[0];
rz(0.068724364) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(1.2083763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2657304) q[0];
sx q[0];
rz(-3.0348572) q[0];
sx q[0];
rz(2.021832) q[0];
rz(-pi) q[1];
rz(-1.5829093) q[2];
sx q[2];
rz(-1.6965869) q[2];
sx q[2];
rz(2.8319401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87885107) q[1];
sx q[1];
rz(-0.95572119) q[1];
sx q[1];
rz(1.4351074) q[1];
rz(-pi) q[2];
rz(2.049202) q[3];
sx q[3];
rz(-2.0432825) q[3];
sx q[3];
rz(0.14379584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(1.6750083) q[2];
rz(2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(-2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0242457) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(0.29719621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8288119) q[0];
sx q[0];
rz(-3.045407) q[0];
sx q[0];
rz(-1.0440774) q[0];
rz(-pi) q[1];
rz(-0.53571312) q[2];
sx q[2];
rz(-2.3539054) q[2];
sx q[2];
rz(2.1802769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6151707) q[1];
sx q[1];
rz(-2.5089658) q[1];
sx q[1];
rz(2.1701865) q[1];
rz(1.1737212) q[3];
sx q[3];
rz(-0.84871549) q[3];
sx q[3];
rz(-1.811036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16195665) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(-2.5276108) q[2];
rz(-2.2654514) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(-0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0456332) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(-0.30763787) q[0];
rz(2.3930507) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(-0.83980733) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9899983) q[0];
sx q[0];
rz(-0.50938207) q[0];
sx q[0];
rz(-0.61181061) q[0];
rz(2.8790881) q[2];
sx q[2];
rz(-0.35048198) q[2];
sx q[2];
rz(-1.0002491) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7786583) q[1];
sx q[1];
rz(-2.0211453) q[1];
sx q[1];
rz(-1.910701) q[1];
rz(-2.967756) q[3];
sx q[3];
rz(-0.94508119) q[3];
sx q[3];
rz(-0.10007773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(-1.9258202) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(-1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76386219) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(2.6960301) q[0];
rz(-2.5207649) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(2.1760118) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27698101) q[0];
sx q[0];
rz(-2.3267713) q[0];
sx q[0];
rz(-0.38397249) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.031164073) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(-0.63809168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.021918745) q[1];
sx q[1];
rz(-1.5389171) q[1];
sx q[1];
rz(1.3171413) q[1];
rz(-pi) q[2];
rz(0.097316381) q[3];
sx q[3];
rz(-1.4574377) q[3];
sx q[3];
rz(1.3219584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.821637) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(2.2122673) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6699566) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(-1.8575645) q[0];
rz(0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(-2.0934385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6402733) q[0];
sx q[0];
rz(-0.87448705) q[0];
sx q[0];
rz(0.93080824) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3088147) q[2];
sx q[2];
rz(-1.1156429) q[2];
sx q[2];
rz(0.27311329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84761274) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(-1.6683679) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1439476) q[3];
sx q[3];
rz(-2.3754658) q[3];
sx q[3];
rz(1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(-0.90083814) q[2];
rz(2.0488996) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(-1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95935217) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(2.545488) q[0];
rz(1.6456564) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398864) q[0];
sx q[0];
rz(-0.6131999) q[0];
sx q[0];
rz(-2.1481811) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2539358) q[2];
sx q[2];
rz(-1.1057901) q[2];
sx q[2];
rz(2.3777865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54904304) q[1];
sx q[1];
rz(-2.2144496) q[1];
sx q[1];
rz(0.93306577) q[1];
rz(-pi) q[2];
rz(0.75002589) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(-2.1961574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(-0.22496741) q[2];
rz(-0.088430017) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46463075) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(2.8616469) q[0];
rz(-1.6784558) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(-2.8889012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8425927) q[0];
sx q[0];
rz(-2.8447066) q[0];
sx q[0];
rz(-1.4070369) q[0];
rz(-1.5964609) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(0.54021013) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.686324) q[1];
sx q[1];
rz(-2.2043921) q[1];
sx q[1];
rz(-0.92647657) q[1];
rz(-pi) q[2];
rz(2.3971862) q[3];
sx q[3];
rz(-2.230847) q[3];
sx q[3];
rz(-2.6403514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(1.6097216) q[2];
rz(1.1931217) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(-2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-3.0367592) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(0.2789467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7689965) q[0];
sx q[0];
rz(-2.0645752) q[0];
sx q[0];
rz(-0.84817024) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2645623) q[2];
sx q[2];
rz(-0.5224723) q[2];
sx q[2];
rz(-0.82447169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5178691) q[1];
sx q[1];
rz(-1.1964799) q[1];
sx q[1];
rz(2.3831297) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99854462) q[3];
sx q[3];
rz(-1.4456985) q[3];
sx q[3];
rz(2.3541114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5489244) q[2];
rz(2.0907949) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46362296) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(-1.2532225) q[0];
rz(0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(-1.1463096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.779594) q[0];
sx q[0];
rz(-2.3906374) q[0];
sx q[0];
rz(0.10303084) q[0];
rz(-pi) q[1];
rz(2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(0.42665542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1100626) q[1];
sx q[1];
rz(-1.4974125) q[1];
sx q[1];
rz(-0.88295464) q[1];
x q[2];
rz(-2.5961612) q[3];
sx q[3];
rz(-0.71119961) q[3];
sx q[3];
rz(1.7033073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(1.0158319) q[2];
rz(0.9097957) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(2.773496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(-0.01424271) q[0];
rz(2.2968538) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(1.3815809) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2505014) q[0];
sx q[0];
rz(-1.4282465) q[0];
sx q[0];
rz(3.13091) q[0];
rz(3.0773452) q[2];
sx q[2];
rz(-2.148743) q[2];
sx q[2];
rz(-0.77677514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8091619) q[1];
sx q[1];
rz(-0.5756439) q[1];
sx q[1];
rz(1.8368506) q[1];
rz(-0.36079447) q[3];
sx q[3];
rz(-0.56305712) q[3];
sx q[3];
rz(2.7891956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(0.98999611) q[2];
rz(-2.2475217) q[3];
sx q[3];
rz(-0.48864135) q[3];
sx q[3];
rz(0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.0902696) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(1.3399711) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(-0.52836616) q[2];
sx q[2];
rz(-2.3522204) q[2];
sx q[2];
rz(1.5632202) q[2];
rz(-0.94974489) q[3];
sx q[3];
rz(-1.2093778) q[3];
sx q[3];
rz(-2.9686684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
