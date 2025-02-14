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
rz(-1.4911574) q[0];
sx q[0];
rz(-2.7661134) q[0];
sx q[0];
rz(-2.7501291) q[0];
rz(0.23671167) q[1];
sx q[1];
rz(4.1799217) q[1];
sx q[1];
rz(7.1290457) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0559056) q[0];
sx q[0];
rz(-2.2330267) q[0];
sx q[0];
rz(-1.7895749) q[0];
x q[1];
rz(-2.691306) q[2];
sx q[2];
rz(-1.992072) q[2];
sx q[2];
rz(2.7305528) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7995092) q[1];
sx q[1];
rz(-1.2984283) q[1];
sx q[1];
rz(1.8567159) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1448946) q[3];
sx q[3];
rz(-2.5509868) q[3];
sx q[3];
rz(-0.27593881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5819431) q[2];
sx q[2];
rz(-1.0407642) q[2];
sx q[2];
rz(2.9259658) q[2];
rz(-0.98254472) q[3];
sx q[3];
rz(-1.6590786) q[3];
sx q[3];
rz(-1.2561579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4172281) q[0];
sx q[0];
rz(-0.094014458) q[0];
sx q[0];
rz(-2.7650058) q[0];
rz(1.4812034) q[1];
sx q[1];
rz(-1.9638289) q[1];
sx q[1];
rz(1.0053763) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621012) q[0];
sx q[0];
rz(-0.32515701) q[0];
sx q[0];
rz(-1.5942911) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3668325) q[2];
sx q[2];
rz(-0.68704359) q[2];
sx q[2];
rz(2.4400638) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1192692) q[1];
sx q[1];
rz(-0.6695153) q[1];
sx q[1];
rz(-0.055107856) q[1];
rz(1.1930258) q[3];
sx q[3];
rz(-2.0888512) q[3];
sx q[3];
rz(-0.6686206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5873831) q[2];
sx q[2];
rz(-2.6770834) q[2];
sx q[2];
rz(2.0655538) q[2];
rz(-0.46766034) q[3];
sx q[3];
rz(-0.83100072) q[3];
sx q[3];
rz(2.9369798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6429546) q[0];
sx q[0];
rz(-2.0641646) q[0];
sx q[0];
rz(2.0306008) q[0];
rz(2.8338762) q[1];
sx q[1];
rz(-1.6650763) q[1];
sx q[1];
rz(1.841338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0205524) q[0];
sx q[0];
rz(-2.5835134) q[0];
sx q[0];
rz(2.9738725) q[0];
rz(0.57794364) q[2];
sx q[2];
rz(-1.5405077) q[2];
sx q[2];
rz(0.89627111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.294178) q[1];
sx q[1];
rz(-1.9929033) q[1];
sx q[1];
rz(0.91011427) q[1];
rz(-2.5505187) q[3];
sx q[3];
rz(-0.58286506) q[3];
sx q[3];
rz(-1.7170563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6230781) q[2];
sx q[2];
rz(-2.0936421) q[2];
sx q[2];
rz(-2.9497362) q[2];
rz(1.1040556) q[3];
sx q[3];
rz(-0.42686978) q[3];
sx q[3];
rz(-1.9328611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0693102) q[0];
sx q[0];
rz(-0.59232124) q[0];
sx q[0];
rz(-2.2533672) q[0];
rz(-0.2725254) q[1];
sx q[1];
rz(-1.2963632) q[1];
sx q[1];
rz(-1.9857508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0391722) q[0];
sx q[0];
rz(-1.2769602) q[0];
sx q[0];
rz(0.0095902082) q[0];
x q[1];
rz(0.49057284) q[2];
sx q[2];
rz(-1.5340048) q[2];
sx q[2];
rz(1.246415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0371746) q[1];
sx q[1];
rz(-1.2942477) q[1];
sx q[1];
rz(2.3642703) q[1];
rz(-pi) q[2];
rz(2.6937595) q[3];
sx q[3];
rz(-1.5671697) q[3];
sx q[3];
rz(1.3161532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3027975) q[2];
sx q[2];
rz(-0.28442997) q[2];
sx q[2];
rz(-0.75308853) q[2];
rz(0.49790844) q[3];
sx q[3];
rz(-1.6585191) q[3];
sx q[3];
rz(2.8411617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4860151) q[0];
sx q[0];
rz(-0.97633728) q[0];
sx q[0];
rz(1.8560386) q[0];
rz(1.9517508) q[1];
sx q[1];
rz(-2.4557476) q[1];
sx q[1];
rz(1.5600342) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6046977) q[0];
sx q[0];
rz(-1.6486537) q[0];
sx q[0];
rz(-1.6998864) q[0];
x q[1];
rz(-1.4194593) q[2];
sx q[2];
rz(-1.8514886) q[2];
sx q[2];
rz(2.7996306) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.469857) q[1];
sx q[1];
rz(-1.7845881) q[1];
sx q[1];
rz(-0.44624568) q[1];
rz(0.26435477) q[3];
sx q[3];
rz(-1.7654382) q[3];
sx q[3];
rz(0.96449145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2994284) q[2];
sx q[2];
rz(-0.93141586) q[2];
sx q[2];
rz(0.13775873) q[2];
rz(2.6970741) q[3];
sx q[3];
rz(-1.3714182) q[3];
sx q[3];
rz(2.3270512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21823068) q[0];
sx q[0];
rz(-2.2683999) q[0];
sx q[0];
rz(2.5905304) q[0];
rz(-2.8463544) q[1];
sx q[1];
rz(-0.72160882) q[1];
sx q[1];
rz(2.3982184) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449042) q[0];
sx q[0];
rz(-1.7040692) q[0];
sx q[0];
rz(0.11479423) q[0];
x q[1];
rz(-0.34140519) q[2];
sx q[2];
rz(-1.7936754) q[2];
sx q[2];
rz(-2.6472732) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8507) q[1];
sx q[1];
rz(-0.85324436) q[1];
sx q[1];
rz(-0.13137551) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4287534) q[3];
sx q[3];
rz(-0.83786115) q[3];
sx q[3];
rz(-0.43213613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4117406) q[2];
sx q[2];
rz(-2.9509632) q[2];
sx q[2];
rz(-0.28212485) q[2];
rz(2.7052346) q[3];
sx q[3];
rz(-1.3151508) q[3];
sx q[3];
rz(2.8859629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2082763) q[0];
sx q[0];
rz(-2.0774807) q[0];
sx q[0];
rz(2.5489885) q[0];
rz(1.1366049) q[1];
sx q[1];
rz(-1.8698144) q[1];
sx q[1];
rz(-1.0677387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48615962) q[0];
sx q[0];
rz(-2.2256333) q[0];
sx q[0];
rz(1.372712) q[0];
x q[1];
rz(0.36892682) q[2];
sx q[2];
rz(-1.4539945) q[2];
sx q[2];
rz(-2.359085) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.021879747) q[1];
sx q[1];
rz(-1.849106) q[1];
sx q[1];
rz(-2.0636419) q[1];
rz(-2.7619902) q[3];
sx q[3];
rz(-2.6983454) q[3];
sx q[3];
rz(-2.8576208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.20256715) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(2.1620046) q[2];
rz(3.0136287) q[3];
sx q[3];
rz(-2.1991859) q[3];
sx q[3];
rz(-1.452182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.99280438) q[0];
sx q[0];
rz(-2.1908741) q[0];
sx q[0];
rz(3.0606781) q[0];
rz(3.0203536) q[1];
sx q[1];
rz(-2.464005) q[1];
sx q[1];
rz(1.1239207) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30643647) q[0];
sx q[0];
rz(-0.84266169) q[0];
sx q[0];
rz(-2.6373498) q[0];
x q[1];
rz(0.40718715) q[2];
sx q[2];
rz(-1.6119119) q[2];
sx q[2];
rz(-2.9585055) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.181001) q[1];
sx q[1];
rz(-0.36534009) q[1];
sx q[1];
rz(2.8613371) q[1];
rz(-1.2865661) q[3];
sx q[3];
rz(-0.54966247) q[3];
sx q[3];
rz(0.46297234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8671882) q[2];
sx q[2];
rz(-0.97325456) q[2];
sx q[2];
rz(-1.6483866) q[2];
rz(-2.6793206) q[3];
sx q[3];
rz(-2.3267764) q[3];
sx q[3];
rz(-1.7260684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129775) q[0];
sx q[0];
rz(-1.7531489) q[0];
sx q[0];
rz(1.9762565) q[0];
rz(-0.041821592) q[1];
sx q[1];
rz(-2.1526497) q[1];
sx q[1];
rz(-1.8080541) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8671252) q[0];
sx q[0];
rz(-2.2750912) q[0];
sx q[0];
rz(2.254266) q[0];
x q[1];
rz(0.012537738) q[2];
sx q[2];
rz(-0.31883966) q[2];
sx q[2];
rz(-0.35976216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4945664) q[1];
sx q[1];
rz(-1.4799282) q[1];
sx q[1];
rz(-0.4161633) q[1];
rz(-1.4708552) q[3];
sx q[3];
rz(-1.458711) q[3];
sx q[3];
rz(-0.76666561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8970783) q[2];
sx q[2];
rz(-1.8316734) q[2];
sx q[2];
rz(-2.6778527) q[2];
rz(-2.5044299) q[3];
sx q[3];
rz(-1.5183247) q[3];
sx q[3];
rz(0.37850982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8265182) q[0];
sx q[0];
rz(-1.0239064) q[0];
sx q[0];
rz(1.6465323) q[0];
rz(2.3118094) q[1];
sx q[1];
rz(-1.2636431) q[1];
sx q[1];
rz(-1.8216546) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6498316) q[0];
sx q[0];
rz(-1.4433089) q[0];
sx q[0];
rz(0.37212917) q[0];
x q[1];
rz(-1.7989203) q[2];
sx q[2];
rz(-1.14883) q[2];
sx q[2];
rz(-2.162148) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5675256) q[1];
sx q[1];
rz(-1.7729771) q[1];
sx q[1];
rz(-0.23781487) q[1];
x q[2];
rz(1.0089031) q[3];
sx q[3];
rz(-1.7911388) q[3];
sx q[3];
rz(-1.6957078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.39184555) q[2];
sx q[2];
rz(-1.4186991) q[2];
sx q[2];
rz(-2.6302272) q[2];
rz(-2.5857207) q[3];
sx q[3];
rz(-1.116773) q[3];
sx q[3];
rz(-0.51978022) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9075254) q[0];
sx q[0];
rz(-3.0341442) q[0];
sx q[0];
rz(-2.7035614) q[0];
rz(3.0829433) q[1];
sx q[1];
rz(-1.6390683) q[1];
sx q[1];
rz(-1.0674089) q[1];
rz(-3.0728523) q[2];
sx q[2];
rz(-2.4394727) q[2];
sx q[2];
rz(-1.7179678) q[2];
rz(1.2642813) q[3];
sx q[3];
rz(-0.96592663) q[3];
sx q[3];
rz(2.929579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
