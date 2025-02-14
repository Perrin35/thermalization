OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.74236953258514) q[0];
sx q[0];
rz(2.72887814243371) q[0];
sx q[0];
rz(5.44698450564548) q[0];
rz(-2.36694812774658) q[1];
sx q[1];
rz(6.34519377549226) q[1];
sx q[1];
rz(11.5580291509549) q[1];
cx q[1],q[0];
rz(1.8275533914566) q[0];
sx q[0];
rz(4.33907476265962) q[0];
sx q[0];
rz(10.9107968568723) q[0];
rz(6.51862239837646) q[2];
sx q[2];
rz(1.94116917450959) q[2];
sx q[2];
rz(10.2067802309911) q[2];
cx q[2],q[1];
rz(-2.85626816749573) q[1];
sx q[1];
rz(3.64597764809663) q[1];
sx q[1];
rz(13.8080944776456) q[1];
rz(3.03622078895569) q[3];
sx q[3];
rz(6.08890763123567) q[3];
sx q[3];
rz(9.05480307935878) q[3];
cx q[3],q[2];
rz(4.44974803924561) q[2];
sx q[2];
rz(4.08512214024598) q[2];
sx q[2];
rz(6.93626401423618) q[2];
rz(3.25609707832336) q[3];
sx q[3];
rz(4.88958183129365) q[3];
sx q[3];
rz(8.36247799395725) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.366505175828934) q[0];
sx q[0];
rz(4.17235687573487) q[0];
sx q[0];
rz(10.7684005260388) q[0];
rz(1.16031813621521) q[1];
sx q[1];
rz(3.55188080866868) q[1];
sx q[1];
rz(11.945787882797) q[1];
cx q[1],q[0];
rz(-1.90898406505585) q[0];
sx q[0];
rz(4.36769345601136) q[0];
sx q[0];
rz(10.7919467449109) q[0];
rz(-0.454994082450867) q[2];
sx q[2];
rz(0.720907600718089) q[2];
sx q[2];
rz(14.2810301542203) q[2];
cx q[2],q[1];
rz(-3.00860404968262) q[1];
sx q[1];
rz(7.06534138520295) q[1];
sx q[1];
rz(10.4706182241361) q[1];
rz(-1.85343897342682) q[3];
sx q[3];
rz(2.16765544016893) q[3];
sx q[3];
rz(9.77386701702281) q[3];
cx q[3],q[2];
rz(1.44299840927124) q[2];
sx q[2];
rz(1.6220761855417) q[2];
sx q[2];
rz(12.7963893175046) q[2];
rz(1.55517852306366) q[3];
sx q[3];
rz(4.2902126630121) q[3];
sx q[3];
rz(9.70046073793575) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.255529373884201) q[0];
sx q[0];
rz(2.75241884787614) q[0];
sx q[0];
rz(8.34390196799442) q[0];
rz(5.63995885848999) q[1];
sx q[1];
rz(3.94094911416108) q[1];
sx q[1];
rz(12.4810554742734) q[1];
cx q[1],q[0];
rz(0.158358305692673) q[0];
sx q[0];
rz(0.933457763987132) q[0];
sx q[0];
rz(7.11781141757175) q[0];
rz(-3.35414242744446) q[2];
sx q[2];
rz(5.27407255967195) q[2];
sx q[2];
rz(14.0121359586637) q[2];
cx q[2],q[1];
rz(1.41156315803528) q[1];
sx q[1];
rz(4.58108249505097) q[1];
sx q[1];
rz(7.68656752108737) q[1];
rz(1.50730764865875) q[3];
sx q[3];
rz(3.76378336747224) q[3];
sx q[3];
rz(8.93183154462978) q[3];
cx q[3],q[2];
rz(-1.20494151115417) q[2];
sx q[2];
rz(3.15115844191099) q[2];
sx q[2];
rz(12.762715792648) q[2];
rz(-0.282285451889038) q[3];
sx q[3];
rz(4.25861969788606) q[3];
sx q[3];
rz(10.4700026273648) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.13337564468384) q[0];
sx q[0];
rz(3.6883928497606) q[0];
sx q[0];
rz(12.9518151044767) q[0];
rz(-4.55499076843262) q[1];
sx q[1];
rz(1.20472672780091) q[1];
sx q[1];
rz(9.17483492790862) q[1];
cx q[1],q[0];
rz(0.378738224506378) q[0];
sx q[0];
rz(2.17024991114671) q[0];
sx q[0];
rz(11.5010030031125) q[0];
rz(0.970863401889801) q[2];
sx q[2];
rz(6.17677346070344) q[2];
sx q[2];
rz(14.6450748205106) q[2];
cx q[2],q[1];
rz(-0.394825965166092) q[1];
sx q[1];
rz(4.54367724259431) q[1];
sx q[1];
rz(4.21626660823032) q[1];
rz(-0.198929741978645) q[3];
sx q[3];
rz(7.07026687462861) q[3];
sx q[3];
rz(13.0732340574185) q[3];
cx q[3],q[2];
rz(-0.461277365684509) q[2];
sx q[2];
rz(1.84674337704713) q[2];
sx q[2];
rz(9.66116846203014) q[2];
rz(1.26050281524658) q[3];
sx q[3];
rz(3.39165571530397) q[3];
sx q[3];
rz(10.1020750164907) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.336097180843353) q[0];
sx q[0];
rz(4.67927983601625) q[0];
sx q[0];
rz(8.93966758846446) q[0];
rz(5.10279607772827) q[1];
sx q[1];
rz(4.68885842164094) q[1];
sx q[1];
rz(7.11968967913791) q[1];
cx q[1],q[0];
rz(2.36199498176575) q[0];
sx q[0];
rz(1.21497240860993) q[0];
sx q[0];
rz(12.229538178436) q[0];
rz(2.73764276504517) q[2];
sx q[2];
rz(4.41821852524812) q[2];
sx q[2];
rz(13.0652119874875) q[2];
cx q[2],q[1];
rz(-0.0665576010942459) q[1];
sx q[1];
rz(4.76039114792878) q[1];
sx q[1];
rz(10.0947115182798) q[1];
rz(-0.0630670711398125) q[3];
sx q[3];
rz(5.11239639123017) q[3];
sx q[3];
rz(7.40850279330417) q[3];
cx q[3],q[2];
rz(5.38466453552246) q[2];
sx q[2];
rz(5.55134764512117) q[2];
sx q[2];
rz(11.7968208551328) q[2];
rz(-0.9292933344841) q[3];
sx q[3];
rz(5.15228167374665) q[3];
sx q[3];
rz(9.19208121895) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.19978511333466) q[0];
sx q[0];
rz(2.65640440781648) q[0];
sx q[0];
rz(10.1936024784963) q[0];
rz(1.35175144672394) q[1];
sx q[1];
rz(3.6146506686979) q[1];
sx q[1];
rz(8.53509083985492) q[1];
cx q[1],q[0];
rz(-7.08055830001831) q[0];
sx q[0];
rz(4.50301328499848) q[0];
sx q[0];
rz(12.1772160291593) q[0];
rz(-1.57214617729187) q[2];
sx q[2];
rz(1.61515155633027) q[2];
sx q[2];
rz(9.00048697590038) q[2];
cx q[2],q[1];
rz(-4.71176385879517) q[1];
sx q[1];
rz(8.0000902732187) q[1];
sx q[1];
rz(14.0902123212735) q[1];
rz(-1.33676850795746) q[3];
sx q[3];
rz(5.82128992875154) q[3];
sx q[3];
rz(14.1062531232755) q[3];
cx q[3],q[2];
rz(6.77348756790161) q[2];
sx q[2];
rz(1.70406654675538) q[2];
sx q[2];
rz(7.83823750018283) q[2];
rz(4.41225385665894) q[3];
sx q[3];
rz(3.56056869228417) q[3];
sx q[3];
rz(9.55680809020206) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.17260575294495) q[0];
sx q[0];
rz(4.27985254128511) q[0];
sx q[0];
rz(7.42726645468875) q[0];
rz(0.830896854400635) q[1];
sx q[1];
rz(1.78327194054658) q[1];
sx q[1];
rz(10.2119087934415) q[1];
cx q[1],q[0];
rz(0.689873039722443) q[0];
sx q[0];
rz(1.70531180699403) q[0];
sx q[0];
rz(10.8629323005597) q[0];
rz(-1.02958393096924) q[2];
sx q[2];
rz(4.10225436289842) q[2];
sx q[2];
rz(11.543914294235) q[2];
cx q[2],q[1];
rz(0.28879189491272) q[1];
sx q[1];
rz(5.18049112160737) q[1];
sx q[1];
rz(13.1018726587216) q[1];
rz(-1.01133227348328) q[3];
sx q[3];
rz(5.04803303082521) q[3];
sx q[3];
rz(10.4440504074018) q[3];
cx q[3],q[2];
rz(2.11810111999512) q[2];
sx q[2];
rz(4.11154160101945) q[2];
sx q[2];
rz(6.2500221490781) q[2];
rz(1.59833323955536) q[3];
sx q[3];
rz(3.85751071770722) q[3];
sx q[3];
rz(4.78110406397983) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.89226603507996) q[0];
sx q[0];
rz(1.58582690556581) q[0];
sx q[0];
rz(8.03167376517459) q[0];
rz(-2.68473672866821) q[1];
sx q[1];
rz(7.42777314980561) q[1];
sx q[1];
rz(8.32419154643222) q[1];
cx q[1],q[0];
rz(-0.83525675535202) q[0];
sx q[0];
rz(3.06557113130624) q[0];
sx q[0];
rz(7.83129796981021) q[0];
rz(3.78750729560852) q[2];
sx q[2];
rz(4.12466249068315) q[2];
sx q[2];
rz(6.93261120318576) q[2];
cx q[2],q[1];
rz(0.173284292221069) q[1];
sx q[1];
rz(4.67846909363801) q[1];
sx q[1];
rz(13.2336044073026) q[1];
rz(-1.24237537384033) q[3];
sx q[3];
rz(5.87202444871003) q[3];
sx q[3];
rz(9.74320564269229) q[3];
cx q[3],q[2];
rz(0.0427723750472069) q[2];
sx q[2];
rz(-0.497511474294118) q[2];
sx q[2];
rz(17.2686953306119) q[2];
rz(0.763573110103607) q[3];
sx q[3];
rz(4.65433517296846) q[3];
sx q[3];
rz(8.38872036933109) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.73246574401855) q[0];
sx q[0];
rz(2.38650730450685) q[0];
sx q[0];
rz(9.71626598238155) q[0];
rz(4.37746334075928) q[1];
sx q[1];
rz(4.33827761014039) q[1];
sx q[1];
rz(6.40608093737766) q[1];
cx q[1],q[0];
rz(-1.93663465976715) q[0];
sx q[0];
rz(2.79111480911309) q[0];
sx q[0];
rz(10.8453598976056) q[0];
rz(0.662465929985046) q[2];
sx q[2];
rz(3.49388328393037) q[2];
sx q[2];
rz(12.9162888288419) q[2];
cx q[2],q[1];
rz(1.11574113368988) q[1];
sx q[1];
rz(-0.234997121495656) q[1];
sx q[1];
rz(12.0742380380551) q[1];
rz(1.57682466506958) q[3];
sx q[3];
rz(1.71541264851625) q[3];
sx q[3];
rz(8.8412542104642) q[3];
cx q[3],q[2];
rz(3.7703320980072) q[2];
sx q[2];
rz(1.74919191201264) q[2];
sx q[2];
rz(11.5028087854306) q[2];
rz(-0.177448034286499) q[3];
sx q[3];
rz(5.32494297822053) q[3];
sx q[3];
rz(11.5480625390927) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.437509566545486) q[0];
sx q[0];
rz(2.6884484906965) q[0];
sx q[0];
rz(4.85211656092807) q[0];
rz(0.600729465484619) q[1];
sx q[1];
rz(6.73298636277253) q[1];
sx q[1];
rz(10.0423151016156) q[1];
cx q[1],q[0];
rz(2.40406322479248) q[0];
sx q[0];
rz(1.69497075875337) q[0];
sx q[0];
rz(9.20750246047183) q[0];
rz(5.54010581970215) q[2];
sx q[2];
rz(1.97361448605592) q[2];
sx q[2];
rz(3.1850361585538) q[2];
cx q[2],q[1];
rz(0.706761956214905) q[1];
sx q[1];
rz(7.32114234765107) q[1];
sx q[1];
rz(4.68803832530185) q[1];
rz(0.997977554798126) q[3];
sx q[3];
rz(1.83169344265992) q[3];
sx q[3];
rz(12.4382700681607) q[3];
cx q[3],q[2];
rz(-4.04286241531372) q[2];
sx q[2];
rz(2.14128855069215) q[2];
sx q[2];
rz(13.5437330961148) q[2];
rz(2.05910015106201) q[3];
sx q[3];
rz(2.26682236989076) q[3];
sx q[3];
rz(9.36925959809824) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.65373432636261) q[0];
sx q[0];
rz(1.0845138152414) q[0];
sx q[0];
rz(10.3739490866582) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3.21161723136902) q[1];
sx q[1];
rz(3.89511272509629) q[1];
sx q[1];
rz(14.1194433927457) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.51753956079483) q[2];
sx q[2];
rz(4.69146135647828) q[2];
sx q[2];
rz(12.0931992292325) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.14433979988098) q[3];
sx q[3];
rz(4.38993087609346) q[3];
sx q[3];
rz(10.1166113376538) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
