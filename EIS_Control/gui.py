from tkinter import Tk, Label, Button



class MainWindow:
    def __init__(self, root):
        self.root = root
        root.title("A simple GUI")
        root.attributes('-topmost', 1)
        root.attributes('-topmost', 0)    

        self.greet_button = Button(root, text="Greet", command=self.greet)
        self.greet_button.grid(row=1)
        
        self.close_button = Button(root, text="Close", command=root.destroy)
        self.close_button.grid(row=1, column=3)
        
        


    def greet(self):
        print("Greetings!")

root = Tk()
my_gui = MainWindow(root)
root.mainloop()