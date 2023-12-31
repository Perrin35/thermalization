OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(-0.36663088) q[0];
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6018574) q[0];
sx q[0];
rz(-2.5849113) q[0];
sx q[0];
rz(2.2599028) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26379649) q[2];
sx q[2];
rz(-0.88636878) q[2];
sx q[2];
rz(2.2823208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75016025) q[1];
sx q[1];
rz(-0.62392985) q[1];
sx q[1];
rz(1.0990012) q[1];
rz(-pi) q[2];
rz(-0.22110181) q[3];
sx q[3];
rz(-3.0924774) q[3];
sx q[3];
rz(-0.55288314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(2.74995) q[2];
rz(-3.0018905) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(-0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249107) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(0.87031594) q[1];
sx q[1];
rz(-1.5753997) q[1];
sx q[1];
rz(1.2044027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28019529) q[0];
sx q[0];
rz(-1.9635927) q[0];
sx q[0];
rz(2.0307721) q[0];
rz(-2.2764678) q[2];
sx q[2];
rz(-1.7127617) q[2];
sx q[2];
rz(2.1324468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(-1.6771392) q[1];
x q[2];
rz(-1.9235839) q[3];
sx q[3];
rz(-2.4376166) q[3];
sx q[3];
rz(-0.17458992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5045972) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(-1.2191999) q[2];
rz(-2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.59042674) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(2.5235126) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(1.9504257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1026099) q[0];
sx q[0];
rz(-1.257574) q[0];
sx q[0];
rz(2.9818176) q[0];
x q[1];
rz(2.4345258) q[2];
sx q[2];
rz(-0.36598772) q[2];
sx q[2];
rz(2.6235839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1515854) q[1];
sx q[1];
rz(-1.3602435) q[1];
sx q[1];
rz(1.3596119) q[1];
rz(-pi) q[2];
rz(-2.5997945) q[3];
sx q[3];
rz(-1.9199315) q[3];
sx q[3];
rz(-0.19402129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96737343) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(1.013914) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(-1.0323662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33092609) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(0.20430918) q[0];
rz(1.3775685) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(-2.2185982) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5343691) q[0];
sx q[0];
rz(-1.2396221) q[0];
sx q[0];
rz(1.8934728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.043945233) q[2];
sx q[2];
rz(-1.5590132) q[2];
sx q[2];
rz(0.95566434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7912485) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(-1.9700325) q[1];
rz(-pi) q[2];
rz(-0.022425671) q[3];
sx q[3];
rz(-2.2750345) q[3];
sx q[3];
rz(-1.1426329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(0.50951177) q[2];
rz(2.4984958) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(-2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5287857) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(2.2221785) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.3580094) q[1];
sx q[1];
rz(1.7808419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6960984) q[0];
sx q[0];
rz(-1.3514263) q[0];
sx q[0];
rz(1.7646837) q[0];
rz(-pi) q[1];
rz(2.4945716) q[2];
sx q[2];
rz(-2.2685452) q[2];
sx q[2];
rz(0.24017142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.67140019) q[1];
sx q[1];
rz(-0.50506401) q[1];
sx q[1];
rz(2.7027674) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1702483) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(-1.5435227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(0.11486593) q[2];
rz(-2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.8615287) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546656) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(1.8967569) q[0];
rz(-0.70760977) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(2.8663666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21101418) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(0.37822322) q[0];
x q[1];
rz(0.49595828) q[2];
sx q[2];
rz(-0.99018103) q[2];
sx q[2];
rz(0.95326391) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6323159) q[1];
sx q[1];
rz(-2.5264611) q[1];
sx q[1];
rz(2.5110911) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1915057) q[3];
sx q[3];
rz(-1.6380777) q[3];
sx q[3];
rz(0.34611191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(-0.7115055) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(-2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28073072) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(-1.3635427) q[0];
rz(-1.4200312) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(-2.4218959) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450711) q[0];
sx q[0];
rz(-2.4356027) q[0];
sx q[0];
rz(2.6939874) q[0];
rz(-pi) q[1];
rz(-2.3709488) q[2];
sx q[2];
rz(-1.0174123) q[2];
sx q[2];
rz(-1.9945952) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81855782) q[1];
sx q[1];
rz(-2.5948988) q[1];
sx q[1];
rz(-0.085303765) q[1];
rz(1.3004488) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(0.28966749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(2.4534295) q[2];
rz(0.041953772) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2280837) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(-1.6465181) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5160617) q[0];
sx q[0];
rz(-2.7357258) q[0];
sx q[0];
rz(2.5748475) q[0];
rz(0.23191339) q[2];
sx q[2];
rz(-0.70966087) q[2];
sx q[2];
rz(2.2556925) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5852768) q[1];
sx q[1];
rz(-1.3116515) q[1];
sx q[1];
rz(-1.6782594) q[1];
rz(-pi) q[2];
rz(-2.5008194) q[3];
sx q[3];
rz(-2.2230004) q[3];
sx q[3];
rz(-0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8994393) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(-3.126826) q[2];
rz(-1.2773369) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167851) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(-2.263608) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(1.3605798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9961062) q[0];
sx q[0];
rz(-1.5331475) q[0];
sx q[0];
rz(2.1426175) q[0];
rz(-2.1979245) q[2];
sx q[2];
rz(-1.6086279) q[2];
sx q[2];
rz(0.19247069) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1189551) q[1];
sx q[1];
rz(-1.9070909) q[1];
sx q[1];
rz(2.789546) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9865773) q[3];
sx q[3];
rz(-1.1559556) q[3];
sx q[3];
rz(-0.14839867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62375162) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(0.8927792) q[2];
rz(-3.1023074) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84353012) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(3.0950586) q[0];
rz(-0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(2.4216901) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7881308) q[0];
sx q[0];
rz(-0.25665584) q[0];
sx q[0];
rz(-1.2528332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5653531) q[2];
sx q[2];
rz(-1.8796225) q[2];
sx q[2];
rz(2.960161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0278922) q[1];
sx q[1];
rz(-0.24735951) q[1];
sx q[1];
rz(-1.5075831) q[1];
rz(-0.7474483) q[3];
sx q[3];
rz(-2.2564853) q[3];
sx q[3];
rz(1.3626584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4204734) q[2];
sx q[2];
rz(-1.7420008) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(-1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-2.9251255) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(0.70710612) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
rz(0.98417102) q[3];
sx q[3];
rz(-1.5979206) q[3];
sx q[3];
rz(0.4978705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
