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
rz(1.7093631029129) q[0];
sx q[0];
rz(4.09945699770982) q[0];
sx q[0];
rz(9.28033312260314) q[0];
rz(0.566751778125763) q[1];
sx q[1];
rz(2.61619386275346) q[1];
sx q[1];
rz(8.44701961278125) q[1];
cx q[1],q[0];
rz(-1.5241961479187) q[0];
sx q[0];
rz(2.86887356837327) q[0];
sx q[0];
rz(8.17592070101901) q[0];
rz(-2.62930035591125) q[2];
sx q[2];
rz(5.70791045029695) q[2];
sx q[2];
rz(12.5775766134183) q[2];
cx q[2],q[1];
rz(-0.493300825357437) q[1];
sx q[1];
rz(6.08365312417085) q[1];
sx q[1];
rz(10.2011726856153) q[1];
rz(-4.89079570770264) q[3];
sx q[3];
rz(3.93690672715242) q[3];
sx q[3];
rz(13.4788236379544) q[3];
cx q[3],q[2];
rz(-0.67291659116745) q[2];
sx q[2];
rz(4.29058578808839) q[2];
sx q[2];
rz(8.49250284432575) q[2];
rz(2.94282341003418) q[3];
sx q[3];
rz(4.25216701825196) q[3];
sx q[3];
rz(8.45941553115054) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.707004904747009) q[0];
sx q[0];
rz(5.37781635125215) q[0];
sx q[0];
rz(12.205243086807) q[0];
rz(1.70656418800354) q[1];
sx q[1];
rz(1.78384927113587) q[1];
sx q[1];
rz(8.60677216052219) q[1];
cx q[1],q[0];
rz(-2.5411102771759) q[0];
sx q[0];
rz(3.62369445164735) q[0];
sx q[0];
rz(7.22752878665134) q[0];
rz(-1.44806110858917) q[2];
sx q[2];
rz(5.081367405253) q[2];
sx q[2];
rz(8.72945002316638) q[2];
cx q[2],q[1];
rz(3.35888385772705) q[1];
sx q[1];
rz(4.27425459225709) q[1];
sx q[1];
rz(8.74864885806247) q[1];
rz(1.51370048522949) q[3];
sx q[3];
rz(1.93400564988191) q[3];
sx q[3];
rz(9.64117666184112) q[3];
cx q[3],q[2];
rz(0.252574175596237) q[2];
sx q[2];
rz(6.73031392891938) q[2];
sx q[2];
rz(11.1457410812299) q[2];
rz(-1.31609261035919) q[3];
sx q[3];
rz(3.89954003890092) q[3];
sx q[3];
rz(9.03654510378047) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.40175676345825) q[0];
sx q[0];
rz(3.63369688590104) q[0];
sx q[0];
rz(11.6954347848813) q[0];
rz(2.82541060447693) q[1];
sx q[1];
rz(6.00162139733369) q[1];
sx q[1];
rz(9.1275616645734) q[1];
cx q[1],q[0];
rz(2.11458277702332) q[0];
sx q[0];
rz(1.50733378727967) q[0];
sx q[0];
rz(10.590435242645) q[0];
rz(0.689156651496887) q[2];
sx q[2];
rz(4.06863233645494) q[2];
sx q[2];
rz(7.56123385428592) q[2];
cx q[2],q[1];
rz(0.199936479330063) q[1];
sx q[1];
rz(4.52572038968141) q[1];
sx q[1];
rz(14.1856956243436) q[1];
rz(-2.23709201812744) q[3];
sx q[3];
rz(0.942164572077342) q[3];
sx q[3];
rz(8.94100911020442) q[3];
cx q[3],q[2];
rz(4.91019105911255) q[2];
sx q[2];
rz(3.96454641421373) q[2];
sx q[2];
rz(6.71078369616672) q[2];
rz(1.18874311447144) q[3];
sx q[3];
rz(3.77154246171052) q[3];
sx q[3];
rz(9.95220926999255) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.4363842010498) q[0];
sx q[0];
rz(1.56428650219972) q[0];
sx q[0];
rz(10.1987025499265) q[0];
rz(0.712908387184143) q[1];
sx q[1];
rz(1.02931013901765) q[1];
sx q[1];
rz(10.1065484046857) q[1];
cx q[1],q[0];
rz(-2.65031576156616) q[0];
sx q[0];
rz(2.80454501708085) q[0];
sx q[0];
rz(6.95935628413364) q[0];
rz(-2.53739619255066) q[2];
sx q[2];
rz(4.89194241364534) q[2];
sx q[2];
rz(8.28407011031314) q[2];
cx q[2],q[1];
rz(0.531827867031097) q[1];
sx q[1];
rz(1.03512874444062) q[1];
sx q[1];
rz(7.34226438998386) q[1];
rz(1.22742366790771) q[3];
sx q[3];
rz(5.05666807492311) q[3];
sx q[3];
rz(7.96231303214237) q[3];
cx q[3],q[2];
rz(4.27437734603882) q[2];
sx q[2];
rz(3.8545499761873) q[2];
sx q[2];
rz(8.46657428740665) q[2];
rz(-2.07516741752625) q[3];
sx q[3];
rz(1.85821476777131) q[3];
sx q[3];
rz(12.1867673158567) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.58500695228577) q[0];
sx q[0];
rz(1.79466942151124) q[0];
sx q[0];
rz(8.86447486876651) q[0];
rz(0.999844610691071) q[1];
sx q[1];
rz(6.07972756226594) q[1];
sx q[1];
rz(11.0467964172284) q[1];
cx q[1],q[0];
rz(-1.04096591472626) q[0];
sx q[0];
rz(4.26914420922334) q[0];
sx q[0];
rz(13.597652411453) q[0];
rz(2.18257212638855) q[2];
sx q[2];
rz(4.24086669285829) q[2];
sx q[2];
rz(9.82264149784251) q[2];
cx q[2],q[1];
rz(2.64743638038635) q[1];
sx q[1];
rz(4.5144146998697) q[1];
sx q[1];
rz(11.7586920022885) q[1];
rz(1.69525110721588) q[3];
sx q[3];
rz(0.961545618372508) q[3];
sx q[3];
rz(6.82746217250034) q[3];
cx q[3],q[2];
rz(1.45242977142334) q[2];
sx q[2];
rz(3.69342956145341) q[2];
sx q[2];
rz(9.20184966026946) q[2];
rz(-3.10685038566589) q[3];
sx q[3];
rz(4.89610997040803) q[3];
sx q[3];
rz(9.4963566198866) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.3361622095108) q[0];
sx q[0];
rz(2.82830888231332) q[0];
sx q[0];
rz(11.4963111638944) q[0];
rz(-1.36097145080566) q[1];
sx q[1];
rz(3.52093428571755) q[1];
sx q[1];
rz(4.99912784098789) q[1];
cx q[1],q[0];
rz(0.4385846555233) q[0];
sx q[0];
rz(4.01395151217515) q[0];
sx q[0];
rz(8.60409650801822) q[0];
rz(0.609833300113678) q[2];
sx q[2];
rz(5.81804648240144) q[2];
sx q[2];
rz(9.05112612842723) q[2];
cx q[2],q[1];
rz(3.32564258575439) q[1];
sx q[1];
rz(4.42473152478273) q[1];
sx q[1];
rz(8.3374290227811) q[1];
rz(0.536216259002686) q[3];
sx q[3];
rz(4.4464660008722) q[3];
sx q[3];
rz(10.797209596626) q[3];
cx q[3],q[2];
rz(0.475089192390442) q[2];
sx q[2];
rz(4.58534267743165) q[2];
sx q[2];
rz(11.9287802934568) q[2];
rz(0.836678862571716) q[3];
sx q[3];
rz(4.24745348294313) q[3];
sx q[3];
rz(8.08072659968539) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.566294312477112) q[0];
sx q[0];
rz(5.85323086579377) q[0];
sx q[0];
rz(9.99231806992694) q[0];
rz(2.71388459205627) q[1];
sx q[1];
rz(4.66908410389955) q[1];
sx q[1];
rz(8.48657039403125) q[1];
cx q[1],q[0];
rz(0.0950312465429306) q[0];
sx q[0];
rz(5.1198195536905) q[0];
sx q[0];
rz(10.6622136592786) q[0];
rz(-0.568455517292023) q[2];
sx q[2];
rz(5.35910669167573) q[2];
sx q[2];
rz(11.3535953521649) q[2];
cx q[2],q[1];
rz(2.77925562858582) q[1];
sx q[1];
rz(4.35650232632691) q[1];
sx q[1];
rz(8.73300877808734) q[1];
rz(-0.863578021526337) q[3];
sx q[3];
rz(1.33613172371919) q[3];
sx q[3];
rz(9.59160186945602) q[3];
cx q[3],q[2];
rz(-2.26188349723816) q[2];
sx q[2];
rz(1.96770134766633) q[2];
sx q[2];
rz(11.1803236961286) q[2];
rz(-1.32276797294617) q[3];
sx q[3];
rz(5.13338461716706) q[3];
sx q[3];
rz(9.39574589057966) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.267385274171829) q[0];
sx q[0];
rz(5.94978681405122) q[0];
sx q[0];
rz(7.71706483363315) q[0];
rz(-1.86776149272919) q[1];
sx q[1];
rz(5.14719096024568) q[1];
sx q[1];
rz(10.256045138828) q[1];
cx q[1],q[0];
rz(-0.481526374816895) q[0];
sx q[0];
rz(2.71140814025933) q[0];
sx q[0];
rz(10.6363702773969) q[0];
rz(-1.44193303585052) q[2];
sx q[2];
rz(5.34336772759492) q[2];
sx q[2];
rz(13.8171295881192) q[2];
cx q[2],q[1];
rz(-0.238493278622627) q[1];
sx q[1];
rz(4.22330990632112) q[1];
sx q[1];
rz(10.210379934303) q[1];
rz(1.34042143821716) q[3];
sx q[3];
rz(4.92996195157106) q[3];
sx q[3];
rz(14.7754173040311) q[3];
cx q[3],q[2];
rz(-1.56358087062836) q[2];
sx q[2];
rz(4.51439491112763) q[2];
sx q[2];
rz(8.46041516064807) q[2];
rz(1.97808051109314) q[3];
sx q[3];
rz(2.53012988169725) q[3];
sx q[3];
rz(6.94211504458591) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.373738557100296) q[0];
sx q[0];
rz(2.40661940177018) q[0];
sx q[0];
rz(10.4021541237752) q[0];
rz(4.89661550521851) q[1];
sx q[1];
rz(4.9770474751764) q[1];
sx q[1];
rz(8.31904063224002) q[1];
cx q[1],q[0];
rz(-0.104066364467144) q[0];
sx q[0];
rz(5.82129731972749) q[0];
sx q[0];
rz(8.11831126212283) q[0];
rz(-2.47491669654846) q[2];
sx q[2];
rz(1.65340009530122) q[2];
sx q[2];
rz(8.65434280633136) q[2];
cx q[2],q[1];
rz(2.21782064437866) q[1];
sx q[1];
rz(5.98491540749604) q[1];
sx q[1];
rz(8.74824432133838) q[1];
rz(-3.42977476119995) q[3];
sx q[3];
rz(1.90872207482392) q[3];
sx q[3];
rz(10.9334918022077) q[3];
cx q[3],q[2];
rz(-4.67487335205078) q[2];
sx q[2];
rz(3.52645108302171) q[2];
sx q[2];
rz(15.8577604055326) q[2];
rz(1.37305653095245) q[3];
sx q[3];
rz(4.87721231778199) q[3];
sx q[3];
rz(11.244915342323) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.493664294481277) q[0];
sx q[0];
rz(2.62390080292756) q[0];
sx q[0];
rz(9.55202216505214) q[0];
rz(1.48089015483856) q[1];
sx q[1];
rz(2.77911853988702) q[1];
sx q[1];
rz(9.25697875618144) q[1];
cx q[1],q[0];
rz(0.0848993360996246) q[0];
sx q[0];
rz(3.17741707538302) q[0];
sx q[0];
rz(10.4302475213925) q[0];
rz(2.17895364761353) q[2];
sx q[2];
rz(4.58505013783509) q[2];
sx q[2];
rz(4.73940846919223) q[2];
cx q[2],q[1];
rz(-0.643931090831757) q[1];
sx q[1];
rz(1.57764581044252) q[1];
sx q[1];
rz(10.4071222305219) q[1];
rz(0.118129417300224) q[3];
sx q[3];
rz(7.97644868691499) q[3];
sx q[3];
rz(11.5694484472196) q[3];
cx q[3],q[2];
rz(-0.0267135128378868) q[2];
sx q[2];
rz(0.938888700800486) q[2];
sx q[2];
rz(10.1859236717145) q[2];
rz(0.0900276973843575) q[3];
sx q[3];
rz(4.14476028283174) q[3];
sx q[3];
rz(13.5169110059659) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.889720141887665) q[0];
sx q[0];
rz(4.4176471551233) q[0];
sx q[0];
rz(8.99355811475917) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-1.35137248039246) q[1];
sx q[1];
rz(8.62446466286714) q[1];
sx q[1];
rz(6.52284905909702) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(5.12816572189331) q[2];
sx q[2];
rz(4.41094938118989) q[2];
sx q[2];
rz(10.1515775680463) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(2.47587490081787) q[3];
sx q[3];
rz(5.39052382309968) q[3];
sx q[3];
rz(11.8561532258908) q[3];
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
