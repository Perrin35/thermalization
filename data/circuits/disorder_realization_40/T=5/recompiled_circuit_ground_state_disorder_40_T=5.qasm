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
rz(0.523878276348114) q[0];
sx q[0];
rz(3.73593059380586) q[0];
sx q[0];
rz(9.11596796511813) q[0];
rz(0.297696948051453) q[1];
sx q[1];
rz(4.39900639851625) q[1];
sx q[1];
rz(10.0367513656537) q[1];
cx q[1],q[0];
rz(0.785917222499847) q[0];
sx q[0];
rz(4.62265244324739) q[0];
sx q[0];
rz(8.85802129506274) q[0];
rz(-0.50030642747879) q[2];
sx q[2];
rz(3.46659013827378) q[2];
sx q[2];
rz(9.23741955160304) q[2];
cx q[2],q[1];
rz(1.87773156166077) q[1];
sx q[1];
rz(3.50404095848138) q[1];
sx q[1];
rz(9.01919213532611) q[1];
rz(1.6314138174057) q[3];
sx q[3];
rz(2.44818863471086) q[3];
sx q[3];
rz(8.07665321826144) q[3];
cx q[3],q[2];
rz(1.16741263866425) q[2];
sx q[2];
rz(4.31933513482148) q[2];
sx q[2];
rz(9.98976740836307) q[2];
rz(-0.510883390903473) q[3];
sx q[3];
rz(3.38039083977277) q[3];
sx q[3];
rz(11.0390761852185) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.53294718265533) q[0];
sx q[0];
rz(3.42263445456559) q[0];
sx q[0];
rz(10.4205684423368) q[0];
rz(-0.587987244129181) q[1];
sx q[1];
rz(3.58150658209855) q[1];
sx q[1];
rz(11.2469154357831) q[1];
cx q[1],q[0];
rz(1.30367660522461) q[0];
sx q[0];
rz(3.63137296040589) q[0];
sx q[0];
rz(8.9004817366521) q[0];
rz(0.327359884977341) q[2];
sx q[2];
rz(4.64078751404817) q[2];
sx q[2];
rz(9.80313560961887) q[2];
cx q[2],q[1];
rz(-0.381825745105743) q[1];
sx q[1];
rz(2.3593215068155) q[1];
sx q[1];
rz(9.71971592902347) q[1];
rz(1.34186816215515) q[3];
sx q[3];
rz(2.64810380538041) q[3];
sx q[3];
rz(9.67070535420581) q[3];
cx q[3],q[2];
rz(1.0413281917572) q[2];
sx q[2];
rz(4.14558974106843) q[2];
sx q[2];
rz(9.44727121888801) q[2];
rz(0.104475393891335) q[3];
sx q[3];
rz(4.74567356904084) q[3];
sx q[3];
rz(10.2952639818112) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.812588810920715) q[0];
sx q[0];
rz(4.17221585114534) q[0];
sx q[0];
rz(7.78655359744235) q[0];
rz(0.0809878557920456) q[1];
sx q[1];
rz(3.80051282246644) q[1];
sx q[1];
rz(10.5948829412381) q[1];
cx q[1],q[0];
rz(0.884298205375671) q[0];
sx q[0];
rz(3.74259785016114) q[0];
sx q[0];
rz(10.7084322929303) q[0];
rz(-0.198624923825264) q[2];
sx q[2];
rz(4.07747963269288) q[2];
sx q[2];
rz(9.60543987750217) q[2];
cx q[2],q[1];
rz(1.24182105064392) q[1];
sx q[1];
rz(4.68373564084107) q[1];
sx q[1];
rz(9.92493436335727) q[1];
rz(0.430960983037949) q[3];
sx q[3];
rz(3.69069746335084) q[3];
sx q[3];
rz(10.1680454969327) q[3];
cx q[3],q[2];
rz(0.662042319774628) q[2];
sx q[2];
rz(2.45526710351045) q[2];
sx q[2];
rz(9.38157112001582) q[2];
rz(0.197338119149208) q[3];
sx q[3];
rz(4.17295268376405) q[3];
sx q[3];
rz(8.99243996142551) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.25468051433563) q[0];
sx q[0];
rz(4.22748282750184) q[0];
sx q[0];
rz(9.09928072094127) q[0];
rz(-0.855728507041931) q[1];
sx q[1];
rz(3.55653909047181) q[1];
sx q[1];
rz(11.4802260160367) q[1];
cx q[1],q[0];
rz(0.0864214897155762) q[0];
sx q[0];
rz(4.22277429898316) q[0];
sx q[0];
rz(9.57069798409148) q[0];
rz(0.735771000385284) q[2];
sx q[2];
rz(5.30895248253877) q[2];
sx q[2];
rz(10.6707876682202) q[2];
cx q[2],q[1];
rz(1.1126743555069) q[1];
sx q[1];
rz(3.68753084738786) q[1];
sx q[1];
rz(9.10408068298503) q[1];
rz(0.143020004034042) q[3];
sx q[3];
rz(4.62650230725343) q[3];
sx q[3];
rz(10.110940194122) q[3];
cx q[3],q[2];
rz(0.934227406978607) q[2];
sx q[2];
rz(4.11725244124467) q[2];
sx q[2];
rz(9.60864647328063) q[2];
rz(1.3217226266861) q[3];
sx q[3];
rz(3.84282121260697) q[3];
sx q[3];
rz(9.29177138804599) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.701652705669403) q[0];
sx q[0];
rz(3.44942480524118) q[0];
sx q[0];
rz(10.3616958022039) q[0];
rz(-0.914947152137756) q[1];
sx q[1];
rz(3.99434688885743) q[1];
sx q[1];
rz(11.1066666603009) q[1];
cx q[1],q[0];
rz(1.10323143005371) q[0];
sx q[0];
rz(3.84425929387147) q[0];
sx q[0];
rz(9.28102164565727) q[0];
rz(0.979165971279144) q[2];
sx q[2];
rz(4.17755595048005) q[2];
sx q[2];
rz(10.0936354756276) q[2];
cx q[2],q[1];
rz(0.256816953420639) q[1];
sx q[1];
rz(3.90172806580598) q[1];
sx q[1];
rz(9.2831125318925) q[1];
rz(0.263033747673035) q[3];
sx q[3];
rz(4.33040741284425) q[3];
sx q[3];
rz(9.81895292400523) q[3];
cx q[3],q[2];
rz(1.09782838821411) q[2];
sx q[2];
rz(4.94731846650178) q[2];
sx q[2];
rz(9.28077802657291) q[2];
rz(1.06756591796875) q[3];
sx q[3];
rz(3.48552867968614) q[3];
sx q[3];
rz(10.116255259506) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.30957925319672) q[0];
sx q[0];
rz(3.94895777304704) q[0];
sx q[0];
rz(10.1931366085927) q[0];
rz(0.857530295848846) q[1];
sx q[1];
rz(4.18808856804902) q[1];
sx q[1];
rz(9.97842291592761) q[1];
cx q[1],q[0];
rz(0.346174329519272) q[0];
sx q[0];
rz(3.66859236558015) q[0];
sx q[0];
rz(9.94625923632785) q[0];
rz(1.34809708595276) q[2];
sx q[2];
rz(3.95527550776536) q[2];
sx q[2];
rz(11.1751156806867) q[2];
cx q[2],q[1];
rz(0.0971156135201454) q[1];
sx q[1];
rz(3.85210117896134) q[1];
sx q[1];
rz(10.1039566159169) q[1];
rz(0.849253594875336) q[3];
sx q[3];
rz(4.24321308930451) q[3];
sx q[3];
rz(9.75164864062473) q[3];
cx q[3],q[2];
rz(1.38626098632813) q[2];
sx q[2];
rz(3.57219720085198) q[2];
sx q[2];
rz(9.87058100699588) q[2];
rz(1.24659943580627) q[3];
sx q[3];
rz(4.91368285019929) q[3];
sx q[3];
rz(11.7163231134336) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0691394731402397) q[0];
sx q[0];
rz(4.10156884987886) q[0];
sx q[0];
rz(9.564849010102) q[0];
rz(0.927992939949036) q[1];
sx q[1];
rz(4.47840658028657) q[1];
sx q[1];
rz(7.97472593783542) q[1];
cx q[1],q[0];
rz(1.01940131187439) q[0];
sx q[0];
rz(3.33572775323922) q[0];
sx q[0];
rz(10.0277080297391) q[0];
rz(1.33947396278381) q[2];
sx q[2];
rz(4.33910325367982) q[2];
sx q[2];
rz(10.5581456184308) q[2];
cx q[2],q[1];
rz(0.191810250282288) q[1];
sx q[1];
rz(5.49666109879548) q[1];
sx q[1];
rz(10.3648164629857) q[1];
rz(2.27367281913757) q[3];
sx q[3];
rz(3.37787656684453) q[3];
sx q[3];
rz(7.27200195788547) q[3];
cx q[3],q[2];
rz(2.10400867462158) q[2];
sx q[2];
rz(3.31984757085378) q[2];
sx q[2];
rz(9.5086490124385) q[2];
rz(-2.15080785751343) q[3];
sx q[3];
rz(4.33058670361573) q[3];
sx q[3];
rz(11.2783242225568) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.246380761265755) q[0];
sx q[0];
rz(4.44773593743379) q[0];
sx q[0];
rz(9.78164581059619) q[0];
rz(1.69959473609924) q[1];
sx q[1];
rz(3.74368903239305) q[1];
sx q[1];
rz(9.43381416275307) q[1];
cx q[1],q[0];
rz(0.624923169612885) q[0];
sx q[0];
rz(3.46185869176919) q[0];
sx q[0];
rz(10.1952984094541) q[0];
rz(0.386308640241623) q[2];
sx q[2];
rz(4.0496393760019) q[2];
sx q[2];
rz(9.08768556117221) q[2];
cx q[2],q[1];
rz(-1.39560925960541) q[1];
sx q[1];
rz(4.44066754181916) q[1];
sx q[1];
rz(11.0904000759046) q[1];
rz(1.65417015552521) q[3];
sx q[3];
rz(3.5862334390455) q[3];
sx q[3];
rz(8.98793259858295) q[3];
cx q[3],q[2];
rz(1.25701999664307) q[2];
sx q[2];
rz(4.48839917977388) q[2];
sx q[2];
rz(10.130293583862) q[2];
rz(-0.269331187009811) q[3];
sx q[3];
rz(4.21804681618745) q[3];
sx q[3];
rz(10.3942381501119) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.893282890319824) q[0];
sx q[0];
rz(3.94694689114625) q[0];
sx q[0];
rz(11.8568887472074) q[0];
rz(-0.642088830471039) q[1];
sx q[1];
rz(5.84178462822969) q[1];
sx q[1];
rz(8.99940926431819) q[1];
cx q[1],q[0];
rz(0.503737211227417) q[0];
sx q[0];
rz(2.20489725668962) q[0];
sx q[0];
rz(9.96404144763156) q[0];
rz(0.587523221969604) q[2];
sx q[2];
rz(3.81066754658753) q[2];
sx q[2];
rz(10.1437201261441) q[2];
cx q[2],q[1];
rz(0.440084189176559) q[1];
sx q[1];
rz(3.59668123920495) q[1];
sx q[1];
rz(10.7891474723737) q[1];
rz(0.996773600578308) q[3];
sx q[3];
rz(4.45492521126802) q[3];
sx q[3];
rz(10.0727961420934) q[3];
cx q[3],q[2];
rz(1.87132585048676) q[2];
sx q[2];
rz(2.96238029201562) q[2];
sx q[2];
rz(8.94598618745013) q[2];
rz(0.350179612636566) q[3];
sx q[3];
rz(4.31944993336732) q[3];
sx q[3];
rz(9.85174060463115) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.66788524389267) q[0];
sx q[0];
rz(3.43852132757241) q[0];
sx q[0];
rz(10.2621254682462) q[0];
rz(0.266520231962204) q[1];
sx q[1];
rz(4.46146038373048) q[1];
sx q[1];
rz(10.7822098493497) q[1];
cx q[1],q[0];
rz(-0.116491794586182) q[0];
sx q[0];
rz(3.22815074225003) q[0];
sx q[0];
rz(9.46458309739038) q[0];
rz(0.488008081912994) q[2];
sx q[2];
rz(3.0072524865442) q[2];
sx q[2];
rz(10.0364660382192) q[2];
cx q[2],q[1];
rz(1.38567936420441) q[1];
sx q[1];
rz(2.22306207020814) q[1];
sx q[1];
rz(9.40906254424855) q[1];
rz(1.24562001228333) q[3];
sx q[3];
rz(4.26996246178681) q[3];
sx q[3];
rz(9.99940041302844) q[3];
cx q[3],q[2];
rz(-0.161676555871964) q[2];
sx q[2];
rz(3.67566797335679) q[2];
sx q[2];
rz(9.86294174789592) q[2];
rz(0.915803730487823) q[3];
sx q[3];
rz(4.66514244874055) q[3];
sx q[3];
rz(8.62159887551471) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.71003079414368) q[0];
sx q[0];
rz(3.77789679368074) q[0];
sx q[0];
rz(9.99955532550021) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.509978711605072) q[1];
sx q[1];
rz(3.46677005489404) q[1];
sx q[1];
rz(11.5849742650907) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.91958224773407) q[2];
sx q[2];
rz(4.40872302849824) q[2];
sx q[2];
rz(9.36318461447164) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.360352367162704) q[3];
sx q[3];
rz(3.66805973847444) q[3];
sx q[3];
rz(10.6879353284757) q[3];
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
