OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(2.6497901) q[0];
sx q[0];
rz(9.2368035) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7029019) q[0];
sx q[0];
rz(-0.54649788) q[0];
sx q[0];
rz(1.370907) q[0];
x q[1];
rz(1.8194524) q[2];
sx q[2];
rz(-2.6373632) q[2];
sx q[2];
rz(-0.76489514) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51812664) q[1];
sx q[1];
rz(-1.4033485) q[1];
sx q[1];
rz(-1.8159588) q[1];
rz(-pi) q[2];
rz(-1.7484619) q[3];
sx q[3];
rz(-1.9068309) q[3];
sx q[3];
rz(-0.48776585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1774896) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(2.5906079) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(-1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(2.205251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3204526) q[0];
sx q[0];
rz(-1.3409412) q[0];
sx q[0];
rz(-1.6605404) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78511946) q[2];
sx q[2];
rz(-1.2650507) q[2];
sx q[2];
rz(0.53346764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3766107) q[1];
sx q[1];
rz(-2.2602343) q[1];
sx q[1];
rz(1.4355684) q[1];
x q[2];
rz(-1.5283269) q[3];
sx q[3];
rz(-2.6385033) q[3];
sx q[3];
rz(2.344362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(-2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(-0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24580978) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(2.202503) q[0];
rz(2.242873) q[1];
sx q[1];
rz(-2.6627314) q[1];
sx q[1];
rz(-2.5476707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.827841) q[0];
sx q[0];
rz(-1.8694287) q[0];
sx q[0];
rz(-1.9065501) q[0];
rz(1.6410286) q[2];
sx q[2];
rz(-0.68968455) q[2];
sx q[2];
rz(0.94674142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2908823) q[1];
sx q[1];
rz(-1.1516654) q[1];
sx q[1];
rz(-0.23371975) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61641683) q[3];
sx q[3];
rz(-0.55509242) q[3];
sx q[3];
rz(0.68716955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(-1.7017986) q[2];
rz(-2.7539608) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(2.7534527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(-0.50278062) q[0];
rz(-2.373383) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(2.3847413) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3319791) q[0];
sx q[0];
rz(-2.6410714) q[0];
sx q[0];
rz(-0.74762263) q[0];
rz(-pi) q[1];
rz(-3.1286131) q[2];
sx q[2];
rz(-1.105096) q[2];
sx q[2];
rz(-2.4109858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6824324) q[1];
sx q[1];
rz(-1.9285413) q[1];
sx q[1];
rz(0.95174241) q[1];
rz(-0.73369153) q[3];
sx q[3];
rz(-1.9568223) q[3];
sx q[3];
rz(-1.8611849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42671529) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(-1.654401) q[2];
rz(-0.58250827) q[3];
sx q[3];
rz(-1.094386) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(-1.2879397) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(1.0505189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462853) q[0];
sx q[0];
rz(-1.8827569) q[0];
sx q[0];
rz(-2.8739268) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6987223) q[2];
sx q[2];
rz(-2.1211229) q[2];
sx q[2];
rz(1.0843104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1861021) q[1];
sx q[1];
rz(-2.0932066) q[1];
sx q[1];
rz(2.2239457) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1220783) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(0.36908484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(2.5081432) q[2];
rz(-1.9472306) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.441992) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(-0.25318405) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(1.4621428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896478) q[0];
sx q[0];
rz(-1.6501353) q[0];
sx q[0];
rz(0.21110714) q[0];
x q[1];
rz(2.2199549) q[2];
sx q[2];
rz(-1.7643133) q[2];
sx q[2];
rz(-0.18093872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0837005) q[1];
sx q[1];
rz(-1.9533227) q[1];
sx q[1];
rz(-1.7544569) q[1];
x q[2];
rz(1.9147647) q[3];
sx q[3];
rz(-1.868639) q[3];
sx q[3];
rz(2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2016466) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(2.3366826) q[2];
rz(-1.9645875) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(2.0578407) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(2.2139363) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(2.129508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9818078) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(2.3083789) q[0];
rz(-0.32832844) q[2];
sx q[2];
rz(-2.7331181) q[2];
sx q[2];
rz(0.35818737) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.183179) q[1];
sx q[1];
rz(-1.5378386) q[1];
sx q[1];
rz(-0.54380137) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4885694) q[3];
sx q[3];
rz(-1.0136908) q[3];
sx q[3];
rz(1.132387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(1.0144462) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6124509) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(-2.0902436) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(-0.28373757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8003214) q[0];
sx q[0];
rz(-2.8386136) q[0];
sx q[0];
rz(3.0269701) q[0];
rz(2.0869414) q[2];
sx q[2];
rz(-2.7661341) q[2];
sx q[2];
rz(1.5309389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9323862) q[1];
sx q[1];
rz(-1.3750409) q[1];
sx q[1];
rz(-0.50435658) q[1];
x q[2];
rz(-2.3582718) q[3];
sx q[3];
rz(-0.82286994) q[3];
sx q[3];
rz(-0.93572039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6909137) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(-1.4452176) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(3.0723363) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(-1.5690631) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7135895) q[0];
sx q[0];
rz(-1.1134976) q[0];
sx q[0];
rz(0.86488117) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9022129) q[2];
sx q[2];
rz(-0.33472543) q[2];
sx q[2];
rz(-0.3604381) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.442765) q[1];
sx q[1];
rz(-2.8932533) q[1];
sx q[1];
rz(0.34766867) q[1];
x q[2];
rz(0.92215718) q[3];
sx q[3];
rz(-1.7382009) q[3];
sx q[3];
rz(-2.5006014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1853603) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(-0.36758962) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(2.1504413) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(0.64176732) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(2.8737601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56775996) q[0];
sx q[0];
rz(-2.0127675) q[0];
sx q[0];
rz(-0.7451591) q[0];
x q[1];
rz(-1.3875302) q[2];
sx q[2];
rz(-1.4472618) q[2];
sx q[2];
rz(2.0545517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99991998) q[1];
sx q[1];
rz(-1.8598286) q[1];
sx q[1];
rz(2.8029867) q[1];
rz(1.0288826) q[3];
sx q[3];
rz(-1.3241323) q[3];
sx q[3];
rz(-0.32113069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3315167) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(-1.6646741) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.1289566) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(2.3616882) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(2.4895888) q[2];
sx q[2];
rz(-2.3163788) q[2];
sx q[2];
rz(-0.083995081) q[2];
rz(-2.4640502) q[3];
sx q[3];
rz(-1.4031706) q[3];
sx q[3];
rz(-1.8085898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
