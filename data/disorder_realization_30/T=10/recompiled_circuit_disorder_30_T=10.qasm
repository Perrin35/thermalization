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
rz(0.491398721933365) q[0];
sx q[0];
rz(2.87706515391404) q[0];
sx q[0];
rz(9.81920827030345) q[0];
rz(0.00616229418665171) q[1];
sx q[1];
rz(2.8013463636213) q[1];
sx q[1];
rz(10.6247876644055) q[1];
cx q[1],q[0];
rz(0.40987503528595) q[0];
sx q[0];
rz(3.09336758975918) q[0];
sx q[0];
rz(9.30709569751426) q[0];
rz(0.627752363681793) q[2];
sx q[2];
rz(3.67643431027467) q[2];
sx q[2];
rz(9.98856351374789) q[2];
cx q[2],q[1];
rz(0.297722727060318) q[1];
sx q[1];
rz(3.75505355198915) q[1];
sx q[1];
rz(10.1294015407483) q[1];
rz(-1.46035706996918) q[3];
sx q[3];
rz(3.52182549436624) q[3];
sx q[3];
rz(11.0013953208844) q[3];
cx q[3],q[2];
rz(1.83469426631927) q[2];
sx q[2];
rz(4.63516679604585) q[2];
sx q[2];
rz(9.13517472743198) q[2];
rz(0.875371932983398) q[3];
sx q[3];
rz(4.14781251748139) q[3];
sx q[3];
rz(9.36506736501261) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.615254759788513) q[0];
sx q[0];
rz(4.00900653203065) q[0];
sx q[0];
rz(9.76876819729015) q[0];
rz(-0.0843317657709122) q[1];
sx q[1];
rz(3.81098941166932) q[1];
sx q[1];
rz(10.7799470186154) q[1];
cx q[1],q[0];
rz(0.387931346893311) q[0];
sx q[0];
rz(3.75616380770738) q[0];
sx q[0];
rz(9.31207427232667) q[0];
rz(2.212153673172) q[2];
sx q[2];
rz(3.34357711871202) q[2];
sx q[2];
rz(8.86655304431125) q[2];
cx q[2],q[1];
rz(1.72057104110718) q[1];
sx q[1];
rz(3.59542271693284) q[1];
sx q[1];
rz(9.27848326265022) q[1];
rz(0.183524280786514) q[3];
sx q[3];
rz(4.63604703743989) q[3];
sx q[3];
rz(9.63696177899047) q[3];
cx q[3],q[2];
rz(0.295586735010147) q[2];
sx q[2];
rz(3.97494790156419) q[2];
sx q[2];
rz(9.96246144770786) q[2];
rz(0.5028355717659) q[3];
sx q[3];
rz(2.71677950223024) q[3];
sx q[3];
rz(10.4352156877439) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.663533449172974) q[0];
sx q[0];
rz(3.67935863335664) q[0];
sx q[0];
rz(10.0656504988591) q[0];
rz(-0.748698532581329) q[1];
sx q[1];
rz(2.13320318062837) q[1];
sx q[1];
rz(10.4898718357007) q[1];
cx q[1],q[0];
rz(0.0708060637116432) q[0];
sx q[0];
rz(3.50097254117066) q[0];
sx q[0];
rz(9.40681761725947) q[0];
rz(0.67133903503418) q[2];
sx q[2];
rz(3.96002182562883) q[2];
sx q[2];
rz(10.405643439285) q[2];
cx q[2],q[1];
rz(1.99484312534332) q[1];
sx q[1];
rz(3.75929770072038) q[1];
sx q[1];
rz(9.48708103074833) q[1];
rz(0.722556293010712) q[3];
sx q[3];
rz(4.52289238770539) q[3];
sx q[3];
rz(9.46656394600078) q[3];
cx q[3],q[2];
rz(0.377252668142319) q[2];
sx q[2];
rz(4.32099464734132) q[2];
sx q[2];
rz(10.1421561002652) q[2];
rz(0.688503682613373) q[3];
sx q[3];
rz(3.76922932465608) q[3];
sx q[3];
rz(9.72231901287242) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.82729572057724) q[0];
sx q[0];
rz(4.34174135525758) q[0];
sx q[0];
rz(9.67446873187228) q[0];
rz(1.0149792432785) q[1];
sx q[1];
rz(3.43782150943811) q[1];
sx q[1];
rz(9.43589648361459) q[1];
cx q[1],q[0];
rz(0.775155901908875) q[0];
sx q[0];
rz(4.30770603020722) q[0];
sx q[0];
rz(9.91880873440906) q[0];
rz(-0.690946996212006) q[2];
sx q[2];
rz(4.00425419409806) q[2];
sx q[2];
rz(8.60947731732532) q[2];
cx q[2],q[1];
rz(0.309184461832047) q[1];
sx q[1];
rz(3.74476722081239) q[1];
sx q[1];
rz(10.0255970120351) q[1];
rz(-0.192418217658997) q[3];
sx q[3];
rz(1.8835860808664) q[3];
sx q[3];
rz(9.92422938942119) q[3];
cx q[3],q[2];
rz(-0.495423853397369) q[2];
sx q[2];
rz(1.88587060769136) q[2];
sx q[2];
rz(9.14168289899036) q[2];
rz(0.663435339927673) q[3];
sx q[3];
rz(3.73475823004777) q[3];
sx q[3];
rz(8.53669664858981) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.111131422221661) q[0];
sx q[0];
rz(2.90093936224515) q[0];
sx q[0];
rz(9.75660932659313) q[0];
rz(0.49452531337738) q[1];
sx q[1];
rz(4.43943396409089) q[1];
sx q[1];
rz(11.2394630670468) q[1];
cx q[1],q[0];
rz(1.06226122379303) q[0];
sx q[0];
rz(3.90789565642411) q[0];
sx q[0];
rz(10.773204779617) q[0];
rz(2.54636740684509) q[2];
sx q[2];
rz(4.00863728125627) q[2];
sx q[2];
rz(9.45681379958197) q[2];
cx q[2],q[1];
rz(1.2111644744873) q[1];
sx q[1];
rz(4.38554886181886) q[1];
sx q[1];
rz(9.66929822265312) q[1];
rz(-0.60527229309082) q[3];
sx q[3];
rz(15*pi/13) q[3];
sx q[3];
rz(10.5047067165296) q[3];
cx q[3],q[2];
rz(1.14226865768433) q[2];
sx q[2];
rz(3.62355983455712) q[2];
sx q[2];
rz(11.3207100391309) q[2];
rz(1.25495314598083) q[3];
sx q[3];
rz(3.98306301434571) q[3];
sx q[3];
rz(11.7281002759854) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.449288368225098) q[0];
sx q[0];
rz(3.14474665269303) q[0];
sx q[0];
rz(10.1062657594602) q[0];
rz(-0.207551255822182) q[1];
sx q[1];
rz(2.678226800757) q[1];
sx q[1];
rz(10.5405436515729) q[1];
cx q[1],q[0];
rz(0.0782431438565254) q[0];
sx q[0];
rz(5.45924249489839) q[0];
sx q[0];
rz(9.99119809865161) q[0];
rz(1.14075028896332) q[2];
sx q[2];
rz(4.08816078503663) q[2];
sx q[2];
rz(8.7564547419469) q[2];
cx q[2],q[1];
rz(0.115684971213341) q[1];
sx q[1];
rz(3.73029908736283) q[1];
sx q[1];
rz(11.0274252653043) q[1];
rz(1.79557847976685) q[3];
sx q[3];
rz(2.45966789324815) q[3];
sx q[3];
rz(9.34401737003728) q[3];
cx q[3],q[2];
rz(-0.212680220603943) q[2];
sx q[2];
rz(4.43843379815156) q[2];
sx q[2];
rz(9.07045078872844) q[2];
rz(-0.319523304700851) q[3];
sx q[3];
rz(4.29416075547273) q[3];
sx q[3];
rz(9.05289987324878) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.00912439823151) q[0];
sx q[0];
rz(3.37289090652997) q[0];
sx q[0];
rz(10.0991239309232) q[0];
rz(1.11224234104156) q[1];
sx q[1];
rz(3.80609682400758) q[1];
sx q[1];
rz(9.9870992064397) q[1];
cx q[1],q[0];
rz(0.0920985043048859) q[0];
sx q[0];
rz(4.72720912297303) q[0];
sx q[0];
rz(10.0930135011594) q[0];
rz(1.85648393630981) q[2];
sx q[2];
rz(3.4865085204416) q[2];
sx q[2];
rz(9.3554594501774) q[2];
cx q[2],q[1];
rz(-0.795314192771912) q[1];
sx q[1];
rz(3.38430565794046) q[1];
sx q[1];
rz(10.3502707242887) q[1];
rz(0.596625626087189) q[3];
sx q[3];
rz(4.64015606244142) q[3];
sx q[3];
rz(9.83334193228885) q[3];
cx q[3],q[2];
rz(1.10186302661896) q[2];
sx q[2];
rz(4.13545695145661) q[2];
sx q[2];
rz(9.76482451557323) q[2];
rz(0.217510536313057) q[3];
sx q[3];
rz(4.08999225695664) q[3];
sx q[3];
rz(9.74346975087329) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.448892205953598) q[0];
sx q[0];
rz(3.09544388775761) q[0];
sx q[0];
rz(9.02833718656703) q[0];
rz(-0.138928264379501) q[1];
sx q[1];
rz(2.68150922854478) q[1];
sx q[1];
rz(10.9461314439695) q[1];
cx q[1],q[0];
rz(0.28946653008461) q[0];
sx q[0];
rz(2.93311009009416) q[0];
sx q[0];
rz(10.6598304271619) q[0];
rz(0.695176243782043) q[2];
sx q[2];
rz(2.55558285315568) q[2];
sx q[2];
rz(10.17445317506) q[2];
cx q[2],q[1];
rz(0.248904511332512) q[1];
sx q[1];
rz(3.92192593415315) q[1];
sx q[1];
rz(10.910935497276) q[1];
rz(0.140942320227623) q[3];
sx q[3];
rz(3.82575771410997) q[3];
sx q[3];
rz(10.5686679840009) q[3];
cx q[3],q[2];
rz(-0.562693178653717) q[2];
sx q[2];
rz(4.19885543187196) q[2];
sx q[2];
rz(9.71910822986766) q[2];
rz(1.13079047203064) q[3];
sx q[3];
rz(4.52019146283204) q[3];
sx q[3];
rz(10.4204206228177) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.493337363004684) q[0];
sx q[0];
rz(4.05678156216676) q[0];
sx q[0];
rz(9.78411377071544) q[0];
rz(0.946114778518677) q[1];
sx q[1];
rz(2.74555358489091) q[1];
sx q[1];
rz(9.69541249274417) q[1];
cx q[1],q[0];
rz(0.76627379655838) q[0];
sx q[0];
rz(3.09410714929039) q[0];
sx q[0];
rz(10.279651439182) q[0];
rz(2.03300404548645) q[2];
sx q[2];
rz(4.0058024247461) q[2];
sx q[2];
rz(9.1246789753358) q[2];
cx q[2],q[1];
rz(-0.860790193080902) q[1];
sx q[1];
rz(4.75809410412843) q[1];
sx q[1];
rz(9.38026362507745) q[1];
rz(0.861895024776459) q[3];
sx q[3];
rz(4.24556020100648) q[3];
sx q[3];
rz(10.1194137692372) q[3];
cx q[3],q[2];
rz(2.31408071517944) q[2];
sx q[2];
rz(3.37118064065511) q[2];
sx q[2];
rz(8.9693379163663) q[2];
rz(0.811962485313416) q[3];
sx q[3];
rz(4.0235157926851) q[3];
sx q[3];
rz(9.97416106461688) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.510460019111633) q[0];
sx q[0];
rz(4.77484348614747) q[0];
sx q[0];
rz(10.164056813709) q[0];
rz(-0.230706810951233) q[1];
sx q[1];
rz(2.67353269656236) q[1];
sx q[1];
rz(9.91967863439723) q[1];
cx q[1],q[0];
rz(0.659428179264069) q[0];
sx q[0];
rz(3.6292551775747) q[0];
sx q[0];
rz(9.34384821950599) q[0];
rz(1.54986000061035) q[2];
sx q[2];
rz(3.46951320965821) q[2];
sx q[2];
rz(9.02638224362537) q[2];
cx q[2],q[1];
rz(1.43292629718781) q[1];
sx q[1];
rz(2.49279436667497) q[1];
sx q[1];
rz(9.60961627065345) q[1];
rz(0.696775197982788) q[3];
sx q[3];
rz(4.63918701012666) q[3];
sx q[3];
rz(9.11644214986964) q[3];
cx q[3],q[2];
rz(0.133592441678047) q[2];
sx q[2];
rz(4.23308161099488) q[2];
sx q[2];
rz(9.75263125299617) q[2];
rz(-0.192063644528389) q[3];
sx q[3];
rz(3.38933459122712) q[3];
sx q[3];
rz(10.4581771850507) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0138542521744967) q[0];
sx q[0];
rz(3.02851198812062) q[0];
sx q[0];
rz(9.96532604693576) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-3.03233313560486) q[1];
sx q[1];
rz(3.39204585750634) q[1];
sx q[1];
rz(13.2104904413144) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.802351951599121) q[2];
sx q[2];
rz(2.46717438300187) q[2];
sx q[2];
rz(9.7100720167081) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.820538282394409) q[3];
sx q[3];
rz(3.64901343186433) q[3];
sx q[3];
rz(9.97769907712146) q[3];
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
