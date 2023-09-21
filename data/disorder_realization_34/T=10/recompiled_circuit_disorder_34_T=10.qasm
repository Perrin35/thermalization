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
rz(-2.76780867576599) q[0];
sx q[0];
rz(5.84355059464509) q[0];
sx q[0];
rz(6.201865887634) q[0];
rz(-2.49131393432617) q[1];
sx q[1];
rz(4.42500165303285) q[1];
sx q[1];
rz(11.7835736036222) q[1];
cx q[1],q[0];
rz(0.622779428958893) q[0];
sx q[0];
rz(2.88400271733338) q[0];
sx q[0];
rz(6.03887293337985) q[0];
rz(-3.48936319351196) q[2];
sx q[2];
rz(8.99024740059907) q[2];
sx q[2];
rz(3.58632371424838) q[2];
cx q[2],q[1];
rz(0.37384569644928) q[1];
sx q[1];
rz(2.00315526326234) q[1];
sx q[1];
rz(7.99509785174533) q[1];
rz(-0.0893931537866592) q[3];
sx q[3];
rz(8.16860357125337) q[3];
sx q[3];
rz(8.91524789332553) q[3];
cx q[3],q[2];
rz(-2.2471182346344) q[2];
sx q[2];
rz(4.14656713803346) q[2];
sx q[2];
rz(12.4505593538205) q[2];
rz(1.5995020866394) q[3];
sx q[3];
rz(3.23789063294465) q[3];
sx q[3];
rz(8.3714717388074) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.887499213218689) q[0];
sx q[0];
rz(3.69112721283967) q[0];
sx q[0];
rz(9.22944264709159) q[0];
rz(5.90814924240112) q[1];
sx q[1];
rz(4.61762860615785) q[1];
sx q[1];
rz(12.8061880826871) q[1];
cx q[1],q[0];
rz(-1.14479064941406) q[0];
sx q[0];
rz(1.72734335263307) q[0];
sx q[0];
rz(8.10295388697788) q[0];
rz(0.300917983055115) q[2];
sx q[2];
rz(4.51031294663484) q[2];
sx q[2];
rz(11.7213234662931) q[2];
cx q[2],q[1];
rz(2.63101959228516) q[1];
sx q[1];
rz(0.782007129984446) q[1];
sx q[1];
rz(3.21629426478549) q[1];
rz(0.715077936649323) q[3];
sx q[3];
rz(6.30477491219575) q[3];
sx q[3];
rz(11.1716480016629) q[3];
cx q[3],q[2];
rz(3.77886915206909) q[2];
sx q[2];
rz(8.62244764168794) q[2];
sx q[2];
rz(8.09492406844302) q[2];
rz(1.79995334148407) q[3];
sx q[3];
rz(4.78426364262635) q[3];
sx q[3];
rz(8.10007140635654) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.277023702859879) q[0];
sx q[0];
rz(5.03638377984101) q[0];
sx q[0];
rz(6.95134637355014) q[0];
rz(1.65026235580444) q[1];
sx q[1];
rz(-2.44900926749175) q[1];
sx q[1];
rz(11.5004541635434) q[1];
cx q[1],q[0];
rz(-1.07201719284058) q[0];
sx q[0];
rz(7.81305042107637) q[0];
sx q[0];
rz(9.17641810177966) q[0];
rz(8.0831298828125) q[2];
sx q[2];
rz(-0.567295877141408) q[2];
sx q[2];
rz(8.37902543543979) q[2];
cx q[2],q[1];
rz(0.909103393554688) q[1];
sx q[1];
rz(5.63815084298188) q[1];
sx q[1];
rz(8.94563824533626) q[1];
rz(1.45077764987946) q[3];
sx q[3];
rz(8.51171747048432) q[3];
sx q[3];
rz(7.42457244395419) q[3];
cx q[3],q[2];
rz(-1.18759095668793) q[2];
sx q[2];
rz(5.05065122445161) q[2];
sx q[2];
rz(9.392230565093) q[2];
rz(-3.50159072875977) q[3];
sx q[3];
rz(4.26821008523042) q[3];
sx q[3];
rz(7.07475826739475) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.00927591323853) q[0];
sx q[0];
rz(7.83864322503144) q[0];
sx q[0];
rz(15.0012569188993) q[0];
rz(-1.93544054031372) q[1];
sx q[1];
rz(2.80010917981202) q[1];
sx q[1];
rz(10.9115620613019) q[1];
cx q[1],q[0];
rz(1.2457754611969) q[0];
sx q[0];
rz(6.72419730027253) q[0];
sx q[0];
rz(13.4322767019193) q[0];
rz(4.23372411727905) q[2];
sx q[2];
rz(3.80807432730729) q[2];
sx q[2];
rz(10.5291396140973) q[2];
cx q[2],q[1];
rz(2.1176016330719) q[1];
sx q[1];
rz(-1.42114338080352) q[1];
sx q[1];
rz(11.2840547323148) q[1];
rz(-0.744385898113251) q[3];
sx q[3];
rz(1.56581691105897) q[3];
sx q[3];
rz(11.0547342061917) q[3];
cx q[3],q[2];
rz(-1.59657454490662) q[2];
sx q[2];
rz(4.03508880932862) q[2];
sx q[2];
rz(11.5578479528348) q[2];
rz(1.09629833698273) q[3];
sx q[3];
rz(7.51181331475312) q[3];
sx q[3];
rz(12.6963174104612) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.03560447692871) q[0];
sx q[0];
rz(0.851770790415355) q[0];
sx q[0];
rz(14.2476062536161) q[0];
rz(1.58850717544556) q[1];
sx q[1];
rz(7.5190636237436) q[1];
sx q[1];
rz(12.5502965211789) q[1];
cx q[1],q[0];
rz(-0.0947293117642403) q[0];
sx q[0];
rz(4.18612733681733) q[0];
sx q[0];
rz(10.1117099285047) q[0];
rz(-1.57913458347321) q[2];
sx q[2];
rz(7.42328992684419) q[2];
sx q[2];
rz(12.3239712476651) q[2];
cx q[2],q[1];
rz(4.02894353866577) q[1];
sx q[1];
rz(3.87013223965699) q[1];
sx q[1];
rz(6.94689366816684) q[1];
rz(1.208491563797) q[3];
sx q[3];
rz(1.56658962567384) q[3];
sx q[3];
rz(10.8818732261579) q[3];
cx q[3],q[2];
rz(5.17391490936279) q[2];
sx q[2];
rz(-1.02919277350371) q[2];
sx q[2];
rz(11.9444801568906) q[2];
rz(1.09709978103638) q[3];
sx q[3];
rz(5.50647655327851) q[3];
sx q[3];
rz(12.3460631132047) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.9428403377533) q[0];
sx q[0];
rz(3.13829136535572) q[0];
sx q[0];
rz(8.5180109500806) q[0];
rz(2.32689070701599) q[1];
sx q[1];
rz(5.59482398827607) q[1];
sx q[1];
rz(11.3416347265164) q[1];
cx q[1],q[0];
rz(0.463661015033722) q[0];
sx q[0];
rz(1.78607979615266) q[0];
sx q[0];
rz(11.9348294496457) q[0];
rz(-1.36250627040863) q[2];
sx q[2];
rz(8.28809467156465) q[2];
sx q[2];
rz(10.6941673517148) q[2];
cx q[2],q[1];
rz(1.86065566539764) q[1];
sx q[1];
rz(-0.608156291646413) q[1];
sx q[1];
rz(9.52130069433852) q[1];
rz(3.03114628791809) q[3];
sx q[3];
rz(8.25294986565644) q[3];
sx q[3];
rz(12.4058646917264) q[3];
cx q[3],q[2];
rz(-2.14270520210266) q[2];
sx q[2];
rz(2.31105163891847) q[2];
sx q[2];
rz(13.5060686826627) q[2];
rz(3.35488486289978) q[3];
sx q[3];
rz(5.942686470347) q[3];
sx q[3];
rz(4.89280841349765) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.16191029548645) q[0];
sx q[0];
rz(4.10611424048478) q[0];
sx q[0];
rz(8.84440759419605) q[0];
rz(1.0549818277359) q[1];
sx q[1];
rz(7.73610511620576) q[1];
sx q[1];
rz(11.865590786926) q[1];
cx q[1],q[0];
rz(1.12829196453094) q[0];
sx q[0];
rz(4.8246707042032) q[0];
sx q[0];
rz(12.1134197473447) q[0];
rz(1.23013472557068) q[2];
sx q[2];
rz(5.94553843339021) q[2];
sx q[2];
rz(7.46985886096164) q[2];
cx q[2],q[1];
rz(6.91093349456787) q[1];
sx q[1];
rz(5.7106763442331) q[1];
sx q[1];
rz(9.25767306088611) q[1];
rz(-3.23550534248352) q[3];
sx q[3];
rz(3.87450799544389) q[3];
sx q[3];
rz(15.3578705549161) q[3];
cx q[3],q[2];
rz(-7.64782524108887) q[2];
sx q[2];
rz(2.83124551375444) q[2];
sx q[2];
rz(15.7301258802335) q[2];
rz(5.53848314285278) q[3];
sx q[3];
rz(5.16610947449739) q[3];
sx q[3];
rz(5.88255188464328) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-3.92843294143677) q[0];
sx q[0];
rz(5.25959602196748) q[0];
sx q[0];
rz(8.00827047824069) q[0];
rz(-1.37570607662201) q[1];
sx q[1];
rz(1.73886516888673) q[1];
sx q[1];
rz(11.6746585130612) q[1];
cx q[1],q[0];
rz(-2.13368463516235) q[0];
sx q[0];
rz(3.93311020930345) q[0];
sx q[0];
rz(12.3440553903501) q[0];
rz(-2.15095210075378) q[2];
sx q[2];
rz(-1.18594249884551) q[2];
sx q[2];
rz(11.9848303556363) q[2];
cx q[2],q[1];
rz(4.6692214012146) q[1];
sx q[1];
rz(10.3375169356638) q[1];
sx q[1];
rz(6.80993435382053) q[1];
rz(0.659595787525177) q[3];
sx q[3];
rz(4.45907131035859) q[3];
sx q[3];
rz(10.0279649853627) q[3];
cx q[3],q[2];
rz(1.65317583084106) q[2];
sx q[2];
rz(4.33307996590669) q[2];
sx q[2];
rz(13.3506128549497) q[2];
rz(2.63583207130432) q[3];
sx q[3];
rz(3.99411067565019) q[3];
sx q[3];
rz(12.7996072530667) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.66556751728058) q[0];
sx q[0];
rz(1.53717556794221) q[0];
sx q[0];
rz(10.146813607208) q[0];
rz(-0.333239138126373) q[1];
sx q[1];
rz(1.94573405583436) q[1];
sx q[1];
rz(11.201399421684) q[1];
cx q[1],q[0];
rz(5.89922618865967) q[0];
sx q[0];
rz(3.8863561471277) q[0];
sx q[0];
rz(10.9344575166623) q[0];
rz(-0.190170928835869) q[2];
sx q[2];
rz(6.70069995720918) q[2];
sx q[2];
rz(10.9078137636106) q[2];
cx q[2],q[1];
rz(-3.06279277801514) q[1];
sx q[1];
rz(5.05567959149415) q[1];
sx q[1];
rz(13.4126302957456) q[1];
rz(1.77026307582855) q[3];
sx q[3];
rz(5.00390532811219) q[3];
sx q[3];
rz(9.19828253089591) q[3];
cx q[3],q[2];
rz(4.17450571060181) q[2];
sx q[2];
rz(1.76207569439942) q[2];
sx q[2];
rz(1.9101247548978) q[2];
rz(3.17565560340881) q[3];
sx q[3];
rz(1.27681997616822) q[3];
sx q[3];
rz(3.77636334895297) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.0547399520874) q[0];
sx q[0];
rz(3.70767167408998) q[0];
sx q[0];
rz(11.0884573221128) q[0];
rz(1.0832896232605) q[1];
sx q[1];
rz(4.54127672513063) q[1];
sx q[1];
rz(10.3929722070615) q[1];
cx q[1],q[0];
rz(-2.46148419380188) q[0];
sx q[0];
rz(3.66323128541047) q[0];
sx q[0];
rz(7.0135433435361) q[0];
rz(3.82609891891479) q[2];
sx q[2];
rz(-2.52865251700347) q[2];
sx q[2];
rz(14.9979595899503) q[2];
cx q[2],q[1];
rz(0.266075164079666) q[1];
sx q[1];
rz(3.49580952723558) q[1];
sx q[1];
rz(8.28948435782596) q[1];
rz(-6.28248310089111) q[3];
sx q[3];
rz(-0.322142211598806) q[3];
sx q[3];
rz(8.85035232304736) q[3];
cx q[3],q[2];
rz(4.21983909606934) q[2];
sx q[2];
rz(5.55531826813752) q[2];
sx q[2];
rz(9.42346567277938) q[2];
rz(-1.14082682132721) q[3];
sx q[3];
rz(4.4585304578119) q[3];
sx q[3];
rz(11.2237760782163) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.664526104927063) q[0];
sx q[0];
rz(6.87546959717805) q[0];
sx q[0];
rz(10.3546904683034) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.93079423904419) q[1];
sx q[1];
rz(0.432827623682567) q[1];
sx q[1];
rz(7.32027337550327) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(3.59943819046021) q[2];
sx q[2];
rz(8.53504482110078) q[2];
sx q[2];
rz(7.30563304423496) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.85778772830963) q[3];
sx q[3];
rz(7.14828983147676) q[3];
sx q[3];
rz(12.1615416765134) q[3];
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